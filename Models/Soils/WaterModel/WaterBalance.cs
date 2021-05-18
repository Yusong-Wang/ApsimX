namespace Models.WaterModel
{
    using APSIM.Shared.Utilities;
    using Interfaces;
    using Models.Core;
    using Models.Soils.Nutrients;
    using Soils;
    using System;
    using System.Collections.Generic;
    using System.Globalization;
    using System.Linq;
    using Newtonsoft.Json;

    /// <summary>
    /// The SoilWater module is a cascading water balance model that owes much to its precursors in 
    /// CERES (Jones and Kiniry, 1986) and PERFECT(Littleboy et al, 1992). 
    /// The algorithms for redistribution of water throughout the soil profile have been inherited from 
    /// the CERES family of models.
    ///
    /// The water characteristics of the soil are specified in terms of the lower limit (ll15), 
    /// drained upper limit(dul) and saturated(sat) volumetric water contents. Water movement is 
    /// described using separate algorithms for saturated or unsaturated flow. It is notable that 
    /// redistribution of solutes, such as nitrate- and urea-N, is carried out in this module.
    ///
    /// Modifications adopted from PERFECT include:
    /// * the effects of surface residues and crop cover on modifying runoff and reducing potential soil evaporation,
    /// * small rainfall events are lost as first stage evaporation rather than by the slower process of second stage evaporation, and
    /// * specification of the second stage evaporation coefficient(cona) as an input parameter, providing more flexibility for describing differences in long term soil drying due to soil texture and environmental effects.
    ///
    /// The module is interfaced with SurfaceOrganicMatter and crop modules so that simulation of the soil water balance 
    /// responds to change in the status of surface residues and crop cover(via tillage, decomposition and crop growth).
    ///
    /// Enhancements beyond CERES and PERFECT include:
    /// * the specification of swcon for each layer, being the proportion of soil water above dul that drains in one day
    /// * isolation from the code of the coefficients determining diffusivity as a function of soil water
    ///   (used in calculating unsaturated flow).Choice of diffusivity coefficients more appropriate for soil type have been found to improve model performance.
    /// * unsaturated flow is permitted to move water between adjacent soil layers until some nominated gradient in 
    ///   soil water content is achieved, thereby accounting for the effect of gravity on the fully drained soil water profile.
    ///
    /// SoilWater is called by APSIM on a daily basis, and typical of such models, the various processes are calculated consecutively. 
    /// This contrasts with models such as SWIM that solve simultaneously a set of differential equations that describe the flow processes.
    /// </summary>
    [ValidParent(ParentType = typeof(Soil))]
    [ViewName("UserInterface.Views.ProfileView")]
    [PresenterName("UserInterface.Presenters.ProfilePresenter")]
    [Serializable]
    public class WaterBalance : ModelCollectionFromResource, ISoilWater
    {
        /// <summary>Link to the soil properties.</summary>
        [Link]
        private Soil soil = null;
        
        /// <summary>Access the soil physical properties.</summary>
        [Link] 
        private IPhysical soilPhysical = null;

        [Link]
        Sample initial = null;

        [Link]
        private ISummary summary = null;

        /// <summary>Link to the lateral flow model.</summary>
        [Link]
        private LateralFlowModel lateralFlowModel = null;

        /// <summary>Link to the runoff model.</summary>
        [Link(Type = LinkType.Child, ByName = true)]
        private RunoffModel runoffModel = null;

        ///// <summary>Link to the saturated flow model.</summary>
        //[Link]
        //private SaturatedFlowModel saturatedFlow = null;

        ///// <summary>Link to the unsaturated flow model.</summary>
        //[Link]
        //private UnsaturatedFlowModel unsaturatedFlow = null;

        /// <summary>Link to the evaporation model.</summary>
        [Link]
        private EvaporationModel evaporationModel = null;

        /// <summary>Link to the water table model.</summary>
        [Link(Type = LinkType.Child, ByName = true)]
        private WaterTableModel waterTableModel = null;

        [Link(ByName = true)]
        ISolute no3 = null;

        [Link(ByName = true)]
        ISolute nh4 = null;

        [Link(ByName = true)]
        ISolute urea = null;

        [Link(ByName = true, IsOptional = true)]
        ISolute cl = null;

        /// <summary>Irrigation information.</summary>
        [NonSerialized]
        private List<IrrigationApplicationType> irrigations;

        /// <summary>Water content (mm).</summary>
        private double[] waterMM;

        /// <summary>Water content (mm/mm).</summary>
        private double[] waterVolumetric;

        /// <summary>Start date for switch to summer parameters for soil water evaporation (dd-mmm)</summary>
        [Units("dd-mmm")]
        [Caption("Summer date")]
        [Description("Start date for switch to summer parameters for soil water evaporation")]
        public string SummerDate { get; set; } = "1-Nov";

        /// <summary>Cummulative soil water evaporation to reach the end of stage 1 soil water evaporation in summer (a.k.a. U)</summary>
        [Bounds(Lower = 0.0, Upper = 40.0)]
        [Units("mm")]
        [Caption("Summer U")]
        [Description("Cummulative soil water evaporation to reach the end of stage 1 soil water evaporation in summer (a.k.a. U)")]
        public double SummerU { get; set; } = 6;

        /// <summary>Drying coefficient for stage 2 soil water evaporation in summer (a.k.a. ConA)</summary>
        [Bounds(Lower = 0.0, Upper = 10.0)]
        [Caption("Summer ConA")]
        [Description("Drying coefficient for stage 2 soil water evaporation in summer (a.k.a. ConA)")]
        public double SummerCona { get; set; } = 3.5;

        /// <summary>Start date for switch to winter parameters for soil water evaporation (dd-mmm)</summary>
        [Units("dd-mmm")]
        [Caption("Winter date")]
        [Description("Start date for switch to winter parameters for soil water evaporation")]
        public string WinterDate { get; set; } = "1-Apr";

        /// <summary>Cummulative soil water evaporation to reach the end of stage 1 soil water evaporation in winter (a.k.a. U).</summary>
        [Bounds(Lower = 0.0, Upper = 10.0)]
        [Units("mm")]
        [Caption("Winter U")]
        [Description("Cummulative soil water evaporation to reach the end of stage 1 soil water evaporation in winter (a.k.a. U).")]
        public double WinterU { get; set; } = 6;

        /// <summary>Drying coefficient for stage 2 soil water evaporation in winter (a.k.a. ConA)</summary>
        [Bounds(Lower = 0.0, Upper = 10.0)]
        [Caption("Winter ConA")]
        [Description("Drying coefficient for stage 2 soil water evaporation in winter (a.k.a. ConA)")]
        public double WinterCona { get; set; } = 2.5;

        /// <summary>Constant in the soil water diffusivity calculation (mm2/day)</summary>
        [Bounds(Lower = 0.0, Upper = 1000.0)]
        [Units("mm2/day")]
        [Caption("Diffusivity constant")]
        [Description("Constant in the soil water diffusivity calculation")]
        public double DiffusConst { get; set; }

        /// <summary>Effect of soil water storage above the lower limit on soil water diffusivity (/mm)</summary>
        [Bounds(Lower = 0.0, Upper = 100.0)]
        [Units("/mm")]
        [Caption("Diffusivity slope")]
        [Description("Effect of soil water storage above the lower limit on soil water diffusivity")]
        public double DiffusSlope { get; set; }

        /// <summary>Fraction of incoming radiation reflected from bare soil</summary>
        [Bounds(Lower = 0.0, Upper = 1.0)]
        [Caption("Albedo")]
        [Description("Fraction of incoming radiation reflected from bare soil")]
        public double Salb { get; set; }

        /// <summary>Runoff Curve Number (CN) for bare soil with average moisture</summary>
        [Bounds(Lower = 1.0, Upper = 100.0)]
        [Caption("CN bare")]
        [Description("Runoff Curve Number (CN) for bare soil with average moisture")]
        public double CN2Bare { get; set; }

        /// <summary>Gets or sets the cn red.</summary>
        [Description("Max. reduction in curve number due to cover")]
        public double CNRed { get; set; } = 20;


        /// <summary>Gets or sets the cn cov.</summary>
        [Description("Cover for max curve number reduction")]
        public double CNCov { get; set; } = 0.8;

        /// <summary>Basal width of the downslope boundary of the catchment for lateral flow calculations (m).</summary>
        [Bounds(Lower = 0.0, Upper = 1.0e8F)]
        [Units("m")]
        [Caption("Basal width")]
        [Description("Basal width of the downslope boundary of the catchment for lateral flow calculations")]
        public double DischargeWidth { get; set; } = 5;

        /// <summary>Catchment area for later flow calculations (m2).</summary>
        [Bounds(Lower = 0.0, Upper = 1.0e8F)]
        [Units("m2")]
        [Caption("Catchment")]
        [Description("Catchment area for lateral flow calculations")]
        public double CatchmentArea { get; set; } = 10;

        /// <summary>Depth strings. Wrapper around Thickness.</summary>
        [JsonIgnore]
        [Description("Depth")]
        [Units("cm")]
        public string[] Depth
        {
            get
            {
                return SoilUtilities.ToDepthStrings(Thickness);
            }
            set
            {
                Thickness = SoilUtilities.ToThickness(value);
            }
        }

        /// <summary>Soil layer thickness for each layer (mm).</summary>
        [Units("mm")]
        [Description("Soil layer thickness for each layer")]
        public double[] Thickness { get; set; }

        /// <summary>Amount of water in the soil (mm).</summary>
        [JsonIgnore]
        public double[] Water 
        { 
            get { return waterMM; } 
            set 
            { 
                waterMM = value; 
                waterVolumetric = MathUtilities.Divide(value, soilPhysical.Thickness); 
            } 
        }

        /// <summary>Amount of water in the soil (mm/mm).</summary>
        [JsonIgnore]
        public double[] SW
        {
            get { return waterVolumetric; }
            set
            {
                waterVolumetric = value;
                waterMM = MathUtilities.Multiply(value, soilPhysical.Thickness);
            }
        }

        /// <summary>Runon (mm).</summary>
        [JsonIgnore]
        public double Runon { get; set; }

        /// <summary>The efficiency (0-1) that solutes move down with water.</summary>
        [JsonIgnore]
        public double[] SoluteFluxEfficiency { get; set; }

        /// <summary>The efficiency (0-1) that solutes move up with water.</summary>
        [JsonIgnore]
        public double[] SoluteFlowEfficiency { get; set; }

        /// <summary> This is set by Microclimate and is rainfall less that intercepted by the canopy and residue components </summary>
        [JsonIgnore]
        public double PotentialInfiltration { get; set; }

        // --- Outputs -------------------------------------------------------------------

        /// <summary>Lateral flow (mm).</summary>
        [JsonIgnore]
        public double[] LateralFlow { get { return lateralFlowModel.OutFlow; } }

        /// <summary>Amount of water moving laterally out of the profile (mm)</summary>
        [JsonIgnore]
        public double[] LateralOutflow { get { return LateralFlow; } }

        /// <summary>Runoff (mm).</summary>
        [JsonIgnore]
        public double Runoff { get; private set; }

        /// <summary>Infiltration (mm).</summary>
        [JsonIgnore]
        public double Infiltration { get; private set; }

        /// <summary>Drainage (mm).</summary>
        [JsonIgnore]
        public double Drainage { get { if (Flux == null) return 0; else return Flux[Flux.Length - 1]; } }

        /// <summary>Evaporation (mm).</summary>
        [JsonIgnore]
        public double Evaporation { get { return evaporationModel.Es; } }

        /// <summary>Water table.</summary>
        [JsonIgnore]
        public double WaterTable { get { return waterTableModel.Depth; } set { waterTableModel.Set(value); } }

        /// <summary>Flux. Water moving down (mm).</summary>
        [JsonIgnore]
        public double[] Flux { get; private set; }

        /// <summary>Flow. Water moving up (mm).</summary>
        [JsonIgnore]
        public double[] Flow { get; private set; }

        /// <summary>Gets todays potential runoff (mm).</summary>
        [JsonIgnore]
        public double PotentialRunoff
        {
            get
            {
                double waterForRunoff = PotentialInfiltration;

                foreach (var irrigation in irrigations)
                {
                    if (irrigation.WillRunoff)
                        waterForRunoff = waterForRunoff + irrigation.Amount;
                }
                return waterForRunoff;
            }
        }

        /// <summary>Provides access to the soil properties.</summary>
        [JsonIgnore]
        public Soil Properties { get { return soil; } }

        ///<summary>Gets soil water content (mm)</summary>
        [JsonIgnore]
        public double[] SWmm { get { return Water; } }

        ///<summary>Gets extractable soil water relative to LL15(mm)</summary>
        [JsonIgnore]
        public double[] ESW { get { return MathUtilities.Subtract(Water, soilPhysical.LL15mm); } }

        ///<summary>Gets potential evaporation from soil surface (mm)</summary>
        [JsonIgnore]
        public double Eos { get { return evaporationModel.Eos; } }

        /// <summary>Gets the actual (realised) soil water evaporation (mm)</summary>
        [JsonIgnore]
        public double Es { get { return evaporationModel.Es; } }

        ///<summary>Time since start of second stage evaporation (days).</summary>
        [JsonIgnore]
        public double T { get { return evaporationModel.t; } }

        /// <summary>Gets potential evapotranspiration of the whole soil-plant system (mm)</summary>
        [JsonIgnore]
        public double Eo { get { return evaporationModel.Eo; } set { evaporationModel.Eo = value; } }

        /// <summary>Fractional amount of water above DUL that can drain under gravity per day.</summary>
        /// <remarks>
        /// Between (SAT and DUL) soil water conductivity constant for each soil layer.
        /// At thicknesses specified in "SoilWater" node of GUI.
        /// Use Soil.SWCON for SWCON in standard thickness
        /// </remarks>
        [Bounds(Lower = 0.0, Upper = 1.0)]
        [Units("/d")]
        [Caption("SWCON")]
        [Description("Fractional amount of water above DUL that can drain under gravity per day (SWCON)")]
        public double[] SWCON { get; set; }

        /// <summary>Lateral saturated hydraulic conductivity (KLAT).</summary>
        /// <remarks>
        /// Lateral flow soil water conductivity constant for each soil layer.
        /// At thicknesses specified in "SoilWater" node of GUI.
        /// Use Soil.KLAT for KLAT in standard thickness
        /// </remarks>
        [Bounds(Lower = 0, Upper = 1.0e3F)]
        [Units("mm/d")]
        [Caption("Klat")]
        [Description("Lateral saturated hydraulic conductivity (KLAT)")]
        public double[] KLAT { get; set; }

        /// <summary>Amount of N leaching as NO3-N from the deepest soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double LeachNO3 { get { if (FlowNO3 == null) return 0; else return FlowNO3.Last(); } }

        /// <summary>Amount of N leaching as NH4-N from the deepest soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double LeachNH4 { get { return 0; } }

        /// <summary>Amount of N leaching as urea-N  from the deepest soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double LeachUrea { get { if (FlowUrea == null) return 0; else return FlowUrea.Last(); } }

        /// <summary>Amount of N leaching as NO3 from each soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double[] FlowNO3 { get; private set; }

        /// <summary>Amount of N leaching as NH4 from each soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double[] FlowNH4 { get; private set; }

        /// <summary>Amount of N leaching as urea from each soil layer (kg /ha)</summary>
        [JsonIgnore]
        public double[] FlowUrea { get; private set; }

        /// <summary> This is set by Microclimate and is rainfall less that intercepted by the canopy and residue components </summary>
        [JsonIgnore]
        public double PrecipitationInterception { get; set; }

        /// <summary>Pond.</summary>
        public double Pond { get { return 0; } }

        /// <summary>Plant available water SW-LL15 (mm/mm).</summary>
        [Units("mm/mm")]
        public double[] PAW
        {
            get
            {
                return APSIM.Shared.APSoil.APSoilUtilities.CalcPAWC(soilPhysical.Thickness,
                                                                  soilPhysical.LL15,
                                                                  SW,
                                                                  null);
            }
        }

        /// <summary>Plant available water SW-LL15 (mm).</summary>
        [Units("mm")]
        public double[] PAWmm
        {
            get
            {
                return MathUtilities.Multiply(PAW, soilPhysical.Thickness);
            }
        }

        // --- Variables for new module --------------------------------------------------

        const double psiad = -1e7;

        /// <summary>Number of soil layers.</summary>
        private int num_layers { get { return soilPhysical.Thickness.Length; } }

        // ISoilHydrology hydraulicModel = new HydraulicModels();
        ISoilHydrology hydraulicModel = new SimpleHydraulicModel();

        // Parameters at quadrature points
        // Consider encapsulate all parameters in a class/struct
        private string QuadratureRule;
        private int num_Qpoints;
        private double[,] QDepth;
        private double[,] QWeight;
        private double[,] QTheta;
        private double[,] QTheta_old;
        private double[,] Qpsi;
        private double[,] QK;
        private double[,] QFlux;

        private double[] NewWater;
        private double[] OldWater;
        private double[] InterfaceFlow;
        private int[] FlowType;

        #region Not in use
        /// <summary>
        /// Gets the hydraulic conductivity at DUL (mm/d)
        /// </summary>
        [Description("Hydraulic conductivity at DUL (mm/d)")]
        [Units("mm/d")]
        [Bounds(Lower = 0.0, Upper = 10.0)]
        public double KDul { get; set; }

        /// <summary>
        /// Gets the matric Potential at DUL (mm)
        /// </summary>
        [Description("Matric Potential at DUL (mm)")]
        [Units("mm")]
        [Bounds(Lower = -1e4, Upper = 0.0)]
        public double PSIDul { get; set; }
        #endregion

        // --- Event handlers ------------------------------------------------------------

        /// <summary>Called when a simulation commences.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The event data.</param>
        [EventSubscribe("Commencing")]
        private void OnStartOfSimulation(object sender, EventArgs e)
        {
            Initialise();

            hydraulicModel.Setup(num_layers, soilPhysical);

            InitCalc();
        }

        /// <summary>Called on start of day.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The event data.</param>
        [EventSubscribe("DoDailyInitialisation")]
        private void OnDoDailyInitialisation(object sender, EventArgs e)
        {
            irrigations.Clear();
            Runon = 0;
        }

        /// <summary>Called when an irrigation occurs.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The event data.</param>
        [EventSubscribe("Irrigated")]
        private void OnIrrigated(object sender, IrrigationApplicationType e)
        {
            irrigations.Add(e);
        }

        /// <summary>Called by CLOCK to let this model do its water movement.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The event data.</param>
        [EventSubscribe("DoSoilWaterMovement")]
        private void OnDoSoilWaterMovement(object sender, EventArgs e)
        {
            // Calculate lateral flow.
            lateralFlowModel.Calculate();
            if (LateralFlow.Length > 0)
                Water = MathUtilities.Subtract(Water, LateralFlow);

            // Calculate runoff.
            Runoff = runoffModel.Value();

            // Calculate infiltration.
            Infiltration = PotentialInfiltration - Runoff;

            // Calculate watermovement
            WaterFlow();

            // Calculate Flux, Runoff, and Water at surface layer. 
            Flux = new double[num_layers];
            for (int layer = 0; layer < num_layers; ++layer)
            {
                Flux[layer] = InterfaceFlow[layer + 1];
            }

            Water[0] = Water[0] + InterfaceFlow[0] + Runon;

            Runoff += Infiltration - InterfaceFlow[0];
            Infiltration = Math.Max(0.0, InterfaceFlow[0]);
                

            // Allow irrigation to infiltrate.
            foreach (var irrigation in irrigations)
            {
                if (irrigation.Amount > 0)
                {
                    int irrigationLayer = SoilUtilities.LayerIndexOfDepth(soilPhysical.Thickness, Convert.ToInt32(irrigation.Depth, CultureInfo.InvariantCulture));
                    if (Convert.ToInt32(irrigation.Depth, CultureInfo.InvariantCulture) != 0)
                        Water[irrigationLayer] += irrigation.Amount;
                    //if (irrigationLayer == 0)
                    //    Infiltration += irrigation.Amount;

                    if (no3 != null)
                        no3.kgha[irrigationLayer] += irrigation.NO3;

                    if (nh4 != null)
                        nh4.kgha[irrigationLayer] += irrigation.NH4;

                    if (cl != null)
                        cl.kgha[irrigationLayer] += irrigation.CL;
                }
            }

            // Saturated flow.
            //Flux = saturatedFlow.Values;

            // Add backed up water to runoff. 
            //Water[0] = Water[0] - saturatedFlow.backedUpSurface;

            // Now reduce the infiltration amount by what backed up.
            //Infiltration = Infiltration - saturatedFlow.backedUpSurface;

            // Turn the proportion of the infiltration that backed up into runoff.
            //Runoff = Runoff + saturatedFlow.backedUpSurface;

            // Should go to pond if one exists.
            //  pond = Math.Min(Runoff, max_pond);
            MoveDown(Water, Flux);
            OldWater = Water;

            double[] no3Values = no3.kgha;
            double[] ureaValues = urea.kgha;

            // Calcualte solute movement down with water.
            double[] no3Down = CalculateSoluteMovementDown(no3Values, Water, Flux, SoluteFluxEfficiency);
            MoveDown(no3Values, no3Down);
            double[] ureaDown = CalculateSoluteMovementDown(ureaValues, Water, Flux, SoluteFluxEfficiency);
            MoveDown(ureaValues, ureaDown);

            // Calculate evaporation and remove from top layer.
            // TODO: temporary remove evaporation for testing.
             double es = evaporationModel.Calculate();
             Water[0] = Water[0] - es;

            // Calculate unsaturated flow of water and apply.
            //Flow = unsaturatedFlow.Values;
            //MoveUp(Water, Flow);

            // Check for errors in water variables.
            //CheckForErrors();

            // Calculate water table depth.
            waterTableModel.Calculate();

            //// Calculate and apply net solute movement.
            //double[] no3Up = CalculateNetSoluteMovement(no3Values, Water, Flow, SoluteFlowEfficiency);
            //MoveUp(no3Values, no3Up);
            //double[] ureaUp = CalculateNetSoluteMovement(ureaValues, Water, Flow, SoluteFlowEfficiency);
            //MoveUp(ureaValues, ureaUp);

            //// Update flow output variables.
            //FlowNO3 = MathUtilities.Subtract(no3Down, no3Up);
            //FlowUrea = MathUtilities.Subtract(ureaDown, ureaUp);

            // Set solute state variables.
            no3.SetKgHa(SoluteSetterType.Soil, no3Values);
            urea.SetKgHa(SoluteSetterType.Soil, ureaValues);

            // Now that we've finished moving water, calculate volumetric water
            waterVolumetric = MathUtilities.Divide(Water, soilPhysical.Thickness);
        }

        /// <summary>Move water down the profile</summary>
        /// <param name="water">The water values</param>
        /// <param name="flux">The amount to move down</param>
        private static void MoveDown(double[] water, double[] flux)
        {
            for (int i = 0; i < water.Length; i++)
            {
                if (i == 0)
                    water[i] = water[i] - flux[i];
                else
                    water[i] = water[i] + flux[i - 1] - flux[i];
            }
        }

        /// <summary>Move water up the profile.</summary>
        /// <param name="water">The water values.</param>
        /// <param name="flow">The amount to move up.</param>
        private static void MoveUp(double[] water, double[] flow)
        {
            for (int i = 0; i < water.Length; i++)
            {
                if (i == 0)
                    water[i] = water[i] + flow[i];
                else
                    water[i] = water[i] + flow[i] - flow[i - 1];
            }
        }

        /// <summary>Calculate the solute movement DOWN based on flux.</summary>
        /// <param name="solute"></param>
        /// <param name="water"></param>
        /// <param name="flux"></param>
        /// <param name="efficiency"></param>
        /// <returns></returns>
        private static double[] CalculateSoluteMovementDown(double[] solute, double[] water, double[] flux, double[] efficiency)
        {
            double[] soluteFlux = new double[solute.Length];
            for (int i = 0; i < solute.Length; i++)
            {
                var soluteInLayer = solute[i];
                if (i > 0)
                    soluteInLayer += soluteFlux[i - 1];

                soluteFlux[i] = soluteInLayer * MathUtilities.Divide(flux[i], water[i] + flux[i], 0) * efficiency[i];
                soluteFlux[i] = MathUtilities.Constrain(soluteFlux[i], 0.0, Math.Max(soluteInLayer, 0));
            }

            return soluteFlux;
        }

        /// <summary>Calculate the solute movement UP and DOWN based on flow.</summary>
        /// <param name="solute"></param>
        /// <param name="water"></param>
        /// <param name="flux"></param>
        /// <param name="efficiency"></param>
        /// <returns></returns>
        private static double[] CalculateNetSoluteMovement(double[] solute, double[] water, double[] flux, double[] efficiency)
        {
            double[] soluteUp = CalculateSoluteMovementUp(solute, water, flux, efficiency);

            double[] remaining = new double[flux.Length];
            remaining[0] = soluteUp[0];
            for (int i = 1; i < solute.Length; i++)
                remaining[i] = soluteUp[i] - soluteUp[i - 1];

            double[] soluteDown = new double[solute.Length];
            for (int i = 0; i < solute.Length; i++)
            {
                if (flux[i] < 0)
                {
                    var positiveFlux = flux[i] * -1;
                    var waterInLayer = water[i] + positiveFlux;
                    var soluteInLayer = solute[i] + remaining[i];
                    if (i > 0)
                    {
                        soluteInLayer += soluteDown[i - 1];
                        waterInLayer += flux[i - 1];
                    }

                    soluteDown[i] = positiveFlux * soluteInLayer / waterInLayer * efficiency[i];
                    soluteDown[i] = MathUtilities.Constrain(soluteDown[i], 0, soluteInLayer);
                }
            }
            return MathUtilities.Subtract(soluteUp, soluteDown);
        }

        /// <summary>Calculate the solute movement UP based on flow.</summary>
        /// <param name="solute"></param>
        /// <param name="water"></param>
        /// <param name="flow"></param>
        /// <param name="efficiency"></param>
        /// <returns></returns>
        private static double[] CalculateSoluteMovementUp(double[] solute, double[] water, double[] flow, double[] efficiency)
        {
            // soluteFlow[i] is the solutes flowing into this layer from the layer below.
            // this is the water moving into this layer * solute concentration. That is,
            // water in this layer * solute in this layer / water in this layer.
            //
            // todo: should this be solute[i + 1] because solute concenctration in the water
            // should actually be the solute concenctration in the water moving into this layer
            // from the layer below.
            // flow[i] is the water coming into a layer from the layer below
            double[] soluteFlow = new double[solute.Length];
            for (int i = solute.Length - 2; i >= 0; i--)
            {
                //if (i == 0)
                //    // soluteFlow[i] = 0;?
                //    soluteFlow[i] = flow[i] * solute[i+1] / (water[i+1] + flow[i]);
                //else if (i < solute.Length-2)
                if (flow[i] <= 0)
                    soluteFlow[i] = 0;
                else
                {
                    var soluteInLayer = solute[i + 1] + soluteFlow[i + 1];
                    soluteFlow[i] = flow[i] * soluteInLayer / (water[i + 1] + flow[i] - flow[i + 1]) * efficiency[i];
                    soluteFlow[i] = MathUtilities.Constrain(soluteFlow[i], 0, soluteInLayer);
                }
            }

            return soluteFlow;
        }

        /// <summary>Checks for soil for errors.</summary>
        private void CheckForErrors()
        {
            const double specific_bd = 2.65;

            double min_sw = 0.0;

            for (int i = 0; i < soilPhysical.Thickness.Length; i++)
            {
                double max_sw = 1.0 - MathUtilities.Divide(soilPhysical.BD[i], specific_bd, 0.0);  // ie. Total Porosity

                if (MathUtilities.IsLessThan(soilPhysical.AirDry[i], min_sw))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4})",
                                               " Air dry lower limit of ",
                                               soilPhysical.AirDry[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is below acceptable value of ",
                                               min_sw));

                if (MathUtilities.IsLessThan(soilPhysical.LL15[i], soilPhysical.AirDry[i]))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4})",
                                               " 15 bar lower limit of ",
                                               soilPhysical.LL15[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is below air dry value of ",
                                               soilPhysical.AirDry[i]));

                if (MathUtilities.IsLessThanOrEqual(soilPhysical.DUL[i], soilPhysical.LL15[i]))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4})",
                                               " drained upper limit of ",
                                               soilPhysical.DUL[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is at or below lower limit of ",
                                               soilPhysical.LL15[i]));

                if (MathUtilities.IsLessThanOrEqual(soilPhysical.SAT[i], soilPhysical.DUL[i]))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4})",
                                               " saturation of ",
                                               soilPhysical.SAT[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is at or below drained upper limit of ",
                                               soilPhysical.DUL[i]));

                if (MathUtilities.IsGreaterThan(soilPhysical.SAT[i], max_sw))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4} {7} {8} {9:G4} {10} {11} {12:G4})",
                                               " saturation of ",
                                               soilPhysical.SAT[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is above acceptable value of ",
                                               max_sw,
                                               "\n",
                                               "You must adjust bulk density (bd) to below ",
                                               (1.0 - soilPhysical.SAT[i]) * specific_bd,
                                               "\n",
                                               "OR saturation (sat) to below ",
                                               max_sw));

                if (MathUtilities.IsGreaterThan(SW[i], soilPhysical.SAT[i]))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4}",
                                               " soil water of ",
                                               SW[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is above saturation of ",
                                               soilPhysical.SAT[i]));

                if (MathUtilities.IsLessThan(SW[i], soilPhysical.AirDry[i]))
                    throw new Exception(String.Format("({0} {1:G4}) {2} {3} {4} {5} {6:G4}",
                                               " soil water of ",
                                               SW[i],
                                               " in layer ",
                                               i,
                                               "\n",
                                               "         is below air-dry value of ",
                                               soilPhysical.AirDry[i]));
            }

        }

        ///<summary>Remove water from the profile</summary>
        public void RemoveWater(double[] amountToRemove)
        {
            Water = MathUtilities.Subtract(Water, amountToRemove);
        }

        /// <summary>Sets the water table.</summary>
        /// <param name="InitialDepth">The initial depth.</param> 
        public void SetWaterTable(double InitialDepth)
        {
            WaterTable = InitialDepth;
        }

        ///<summary>Perform a reset</summary>
        public void Reset()
        {
            summary.WriteMessage(this, "Resetting Soil Water Balance");
            Initialise();
        }

        /// <summary>Initialise the model.</summary>
        private void Initialise()
        {
            FlowNH4 = MathUtilities.CreateArrayOfValues(0.0, Thickness.Length);
            SoluteFlowEfficiency = MathUtilities.CreateArrayOfValues(1.0, Thickness.Length);
            SoluteFluxEfficiency = MathUtilities.CreateArrayOfValues(1.0, Thickness.Length);
            Water = initial.SWmm;
            Runon = 0;
            Runoff = 0;
            PotentialInfiltration = 0;
            Flux = null;
            Flow = null;
            evaporationModel.Initialise();
            irrigations = new List<IrrigationApplicationType>();
        }

        ///<summary>Perform tillage</summary>
        public void Tillage(TillageType Data)
        {
            if ((Data.cn_red <= 0) || (Data.cn_rain <= 0))
            {
                string message = "tillage:- " + Data.Name + " has incorrect values for " + Environment.NewLine +
                    "CN reduction = " + Data.cn_red + Environment.NewLine + "Acc rain     = " + Data.cn_red;
                throw new Exception(message);
            }

            double reduction = MathUtilities.Constrain(Data.cn_red, 0.0, CN2Bare);

            runoffModel.TillageCnCumWater = Data.cn_rain;
            runoffModel.TillageCnRed = reduction;
            runoffModel.CumWaterSinceTillage = 0.0;

            var line = string.Format("Soil tilled. CN reduction = {0}. Cumulative rain = {1}", 
                                     reduction, Data.cn_rain);
            summary.WriteMessage(this, line);
        }

        ///<summary>Perform tillage</summary>
        public void Tillage(string tillageType)
        {
            throw new NotImplementedException();
        }


        ///<summary>Perform initial calculations for hydraulic curves</summary>
        private void InitCalc()
        {
            #region Debug code for hydraulic functions
            // ***Check calculated hydraulic parameters with van Genuchten model.
            //double h = -100.0;

            //do
            //{
            //    double theta;
            //    double K;
            //    string str;
            //    for (int layer = 0; layer < num_layers; ++layer)
            //    {
            //        theta = hydraulicModel.get_theta(layer, h);
            //        K = hydraulicModel.get_K(layer, h);
            //        System.Console.WriteLine("layer: " + layer + ", h: " + h + ", theta: " + theta + ", K: " + K);
            //    }

            //    System.Console.WriteLine("Enter another h or enter a negative number to exit.");
            //    str = System.Console.ReadLine();
            //    h = -Convert.ToDouble(str);
            //} while (h <= 0.0);


            //string str;
            //for (int i = 0; i < 6; ++i)
            //{
            //    double theta;
            //    double K;
            //    //System.Console.WriteLine("Input a h value: ");
            //    //str = System.Console.ReadLine();
            //    //h = Convert.ToDouble(str);
            //    theta = hydraulicModel.get_theta(0, h);
            //    K = hydraulicModel.get_K(0, h);
            //    System.Console.WriteLine("h: " + h + ", theta: " + theta + ", K: " + K);
            //    h -= 200;
            //}
            //double Theta = 0.42;
            //for (int i = 0; i < 6; ++i)
            //{
            //    double hh;
            //    //System.Console.WriteLine("Input a theta value: ");
            //    //str = System.Console.ReadLine();
            //    //Theta = Convert.ToDouble(str);
            //    hh = hydraulicModel.get_h(0, Theta);
            //    System.Console.WriteLine("Theta: " + Theta + ", h: " + hh);
            //    Theta -= 0.02;
            //}
            #endregion

            InterfaceFlow = new double[num_layers + 1];
            FlowType = new int[num_layers + 1];
            NewWater = new double[num_layers];
            OldWater = new double[num_layers];

            OldWater = Water;

            QuadratureRule = "Gaussian";
            num_Qpoints = 3;

            InitQuadrature(QuadratureRule, num_Qpoints);       
        }


        /// <summary>
        /// Setup quadrature points for each layer and initialise them.
        /// </summary>
        private void InitQuadrature(string quadrature_rule = "Gaussian", int num_points = 3)
        {
            double[] y;
            double[] weight;
            y = new double[num_points + 1];
            weight = new double[num_points + 1];

            QDepth = new double[num_layers, num_points + 1];
            QWeight = new double[num_layers, num_points + 1];
            QTheta = new double[num_layers, num_points + 1];
            QTheta_old = new double[num_layers, num_points + 1];
            QK = new double[num_layers, num_points + 1];
            Qpsi = new double[num_layers, num_points + 1];
            QFlux = new double[num_layers, num_Qpoints - 1];

            // Point 0 is the middle point representing the average values of the layer.

            if (quadrature_rule == "Gaussian")
            {
                switch (num_points)
                {
                    case 0:
                        System.Console.WriteLine("Number of quadrature points cannot be 0.");
                        break;
                    case 1:
                        y[1] = 1.0 / 2.0;
                        weight[1] = 1.0;
                        y[0] = 1.0 / 2.0;
                        weight[0] = 1.0;
                        break;
                    case 2:
                        y[0] = 1.0 / 2.0;
                        weight[0] = 1.0;
                        y[1] = (1.0 - 1.0 / Math.Sqrt(3.0)) / 2.0;
                        y[2] = (1.0 + 1.0 / Math.Sqrt(3.0)) / 2.0;
                        weight[1] = 1.0 / 2.0;
                        weight[2] = 1.0 / 2.0;
                        break;
                    case 3:
                        y[0] = 1.0 / 2.0;
                        weight[0] = 1.0;
                        y[1] = (1 - Math.Sqrt(3.0 / 5.0)) / 2.0;
                        y[2] = 1.0 / 2.0;
                        y[3] = (1 + Math.Sqrt(3.0 / 5.0)) / 2.0;
                        weight[1] = 5.0 / 9.0 / 2.0;
                        weight[2] = 8.0 / 9.0 / 2.0;
                        weight[3] = 5.0 / 9.0 / 2.0;
                        break;
                    default:
                        System.Console.WriteLine("Number of quadrature points" + num_points + "not implemented.");
                        break;
                }

                for (int layer = 0; layer < num_layers; ++layer)
                {
                    for (int point = 0; point <= num_points; ++point)
                    {
                        QTheta[layer, point] = SW[layer];
                        QTheta_old[layer, point] = QTheta[layer, point];
                        Qpsi[layer, point] = hydraulicModel.get_h(layer, QTheta[layer, point]);
                        QK[layer, point] = hydraulicModel.get_K(layer, Qpsi[layer, point]);
                        QDepth[layer, point] = soilPhysical.ThicknessCumulative[layer] + (y[point] - 1) * soilPhysical.Thickness[layer];
                        QWeight[layer, point] = weight[point];
                    }
                }
            }
            else
                System.Console.WriteLine("Other quadrature rules not implemented.");
        }

        /// <summary> Calculate water flow at each interface.</summary>
        private void WaterFlow()
        {
            // TODO: include Pond calculation.
            double[] Source = new double[num_layers];
            double[] BackFlow = new double[num_layers];
            double[] ExcessW = new double[num_layers];
            double[] PFlow = new double[num_layers + 1];
            bool[] is_Saturated = new bool[num_layers];
            // int MaxIteration = 100;
            double MaxFlow;
            // double UpFlow;
            // double DownFlow;
            // double dist;
            // double psi_diff;

            for (int layer = 0; layer < num_layers; ++layer)
            {
                NewWater[layer] = Water[layer];
                is_Saturated[layer] = false;
                Source[layer] = 0.0;
                BackFlow[layer] = 0.0;
                // Adjust soil water for each quadrature point if Water was changed outside of this class (root uptake).
                Redistribute(layer);
            }
            
            for (int n_interface = 0; n_interface <= num_layers; ++n_interface)
            {
                FlowType[n_interface] = 0;
                InterfaceFlow[n_interface] = 0.0;
            }

            // Check for sources in the profile, e.g. subsurface irrigation
            foreach (var irrigation in irrigations)
            {
                if (irrigation.Amount > 0)
                {
                    int irrigationLayer = SoilUtilities.LayerIndexOfDepth(soilPhysical.Thickness, Convert.ToInt32(irrigation.Depth, CultureInfo.InvariantCulture));
                    if (Convert.ToInt32(irrigation.Depth, CultureInfo.InvariantCulture) == 0)
                        Infiltration += irrigation.Amount;
                    else
                    {
                        Source[irrigationLayer] += irrigation.Amount;
                    }
                }
            }

            // TODO: compare infiltration with potential evaporation and subtract one from another
            // If net value is evaporation, FlowType[0] = -1; otherwise, FlowType[0] = 1;
            FlowType[0] = 1;

            // TODO: free drainage boundary assumed for the bottom; add other possibilities.

            for (int layer = 0; layer < num_layers; ++layer)
            {
                if (Source[layer] > 0.0)
                {
                    NewWater[layer] += Source[layer];
                    is_Saturated[layer] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer] / soilPhysical.SATmm[layer], 0.9999);
                    ExcessW[layer] = Math.Max(0.0, NewWater[layer] - soilPhysical.SATmm[layer]);
                    Redistribute(layer, Source[layer], 0);
                }
            }

            // Forced flow: surface infiltration, subsurface irrigation (flux type), bottom flux flow.
            #region Forced flow

            if (Infiltration > 0.0)
            {
                MaxFlow = QK[0, 0] * ((Pond - (Qpsi[0, 0]) / QDepth[0, 0] + 1));
                MaxFlow = Math.Max(MaxFlow, 2 * soilPhysical.KS[0]);
                InterfaceFlow[0] = Math.Min(MaxFlow, Infiltration);
                NewWater[0] += InterfaceFlow[0];
                Redistribute(0, InterfaceFlow[0], -1);
                // include the following statements in Redistribute()
                if (ExcessW[0] > 0.0)
                    ExcessW[0] += InterfaceFlow[0];
                else
                {
                    is_Saturated[0] = MathUtilities.IsGreaterThanOrEqual(NewWater[0] / soilPhysical.SATmm[0], 0.9999);
                    ExcessW[0] = Math.Max(0.0, NewWater[0] - soilPhysical.SATmm[0]);
                }
            }

            for (int layer = 0; layer < num_layers - 1; ++layer)
            {
                if (ExcessW[layer] > 0.0)
                {
                    FlowType[layer + 1] = 1;

                    MaxFlow = 2 * soilPhysical.KS[layer + 1];

                    if (ExcessW[layer] > MaxFlow)
                    {
                        FlowType[layer + 1] = 2;
                        InterfaceFlow[layer + 1] = MaxFlow;
                        BackFlow[layer] = ExcessW[layer] - MaxFlow;
                        ExcessW[layer] = BackFlow[layer];
                    }
                    else
                    {
                        InterfaceFlow[layer + 1] = ExcessW[layer];
                        if (ExcessW[layer] > MaxFlow / 2.0)
                            FlowType[layer + 1] = 2;
                    }

                    NewWater[layer] -= InterfaceFlow[layer + 1];
                    NewWater[layer + 1] += InterfaceFlow[layer + 1];
                    if (ExcessW[layer + 1] > 0.0)
                        ExcessW[layer + 1] += InterfaceFlow[layer + 1];
                    else
                    {
                        is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1] / soilPhysical.SATmm[layer + 1], 0.9999);
                        ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                        Redistribute(layer + 1, InterfaceFlow[layer + 1], -1);
                    }

                    // Assume excess water only moves downward when layer of the subsurface irrigation is saturated.
                }
            }

            if (ExcessW[num_layers - 1] > 0.0)
            {
                switch (FlowType[num_layers])
                {
                    case -1:
                        // Upward flux boundary for bottom boundary.
                        // TODO: implement this boundary type.
                        break;

                    case 0:
                        // Free drainage
                        MaxFlow = soilPhysical.KS[num_layers - 1];
                        if (ExcessW[num_layers - 1] > MaxFlow)
                        {
                            FlowType[num_layers] = 2;
                            InterfaceFlow[num_layers] = MaxFlow;
                            BackFlow[num_layers - 1] = ExcessW[num_layers - 1] - MaxFlow;
                            ExcessW[num_layers - 1] = BackFlow[num_layers - 1];
                        }
                        else
                        {
                            InterfaceFlow[num_layers] = ExcessW[num_layers - 1];
                        }
                        NewWater[num_layers - 1] -= InterfaceFlow[num_layers];

                        break;

                    case 1:
                        // Downward flux boundary.
                        break;

                    default:
                        System.Console.WriteLine("This bottom boundary type: " + FlowType[num_layers] + " has not been implemented.");
                        break;
                }
            }

            // Backup flow
            for (int layer = num_layers - 1; layer > 0; --layer)
            {
                if (BackFlow[layer] > 0.0)
                {
                    if (FlowType[layer] > 0)
                    {
                        InterfaceFlow[layer] -= BackFlow[layer];

                        if (InterfaceFlow[layer] < 0.0)
                            FlowType[layer] = -2;
                        else
                            FlowType[layer] = 2;
                    }
                    else
                    {
                        InterfaceFlow[layer] = -BackFlow[layer];
                        FlowType[layer] = -1;
                    }

                    NewWater[layer] -= BackFlow[layer];
                    ExcessW[layer] -= BackFlow[layer];
                    NewWater[layer - 1] += BackFlow[layer];

                    if (BackFlow[layer - 1] > 0.0)
                    {
                        BackFlow[layer - 1] += BackFlow[layer];
                        ExcessW[layer - 1] += BackFlow[layer];
                    }
                    else
                    {
                        Redistribute(layer - 1, -BackFlow[layer], 1);
                        is_Saturated[layer - 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer - 1] / soilPhysical.SATmm[layer - 1], 0.9999);
                        BackFlow[layer - 1] = Math.Max(0.0, NewWater[layer - 1] - soilPhysical.SATmm[layer - 1]);
                        ExcessW[layer - 1] = BackFlow[layer - 1];
                    }
                }
            }

            if (BackFlow[0] > 0.0)
            {
                InterfaceFlow[0] -= BackFlow[0];

                NewWater[0] -= BackFlow[0];
                ExcessW[0] -= BackFlow[0];
            }

            // Excess water should be 0 at the end.
            for (int layer = 0; layer < num_layers; ++layer)
            {
                if (!MathUtilities.FloatsAreEqual(ExcessW[layer], 0.0, 1e-6))
                    System.Console.WriteLine("Excess water is not 0 after forced flow.");
            }

            #endregion


            #region Stepwise drainage

            if (FlowType[num_layers] == 0)
            {
                int num_steps = 5;
                double available;
                double flow;
                double K_initial;
                double theta;
                double[] psi_steps = new double[num_steps];
                double[,] theta_steps = new double[num_steps, num_layers];
                double[,] K_steps = new double[num_steps, num_layers];

                PSIDul = -3400.0;
                KDul = 0.1;

                for (int step = 0; step < num_steps; ++step)
                {
                    psi_steps[step] = PSIDul / Math.Pow(2, num_steps - step - 1);
                    for (int layer = 0; layer < num_layers; ++layer)
                    {
                        theta_steps[step, layer] = hydraulicModel.get_theta(layer, psi_steps[step]);
                        K_steps[step, layer] = hydraulicModel.get_K(layer, psi_steps[step]);
                    }
                }
                for(int step = num_steps - 1; step >=0; --step)
                {
                    for (int layer = 0; layer < num_layers; ++layer)
                    {
                        if (step != 0)
                        {
                            K_steps[step, layer] = (K_steps[step - 1, layer] + K_steps[step, layer]) / 2;
                        }
                        else
                        {
                            K_steps[step, layer] = (soilPhysical.KS[layer] + K_steps[step, layer]) / 2;
                        }
                    }
                }

                for (int step = 0; step < num_steps; ++step)
                {
                    for (int layer = 0; layer < num_layers - 1; ++layer)
                    {
                        if (FlowType[layer + 1] > 2 || FlowType[layer + 1] < 0)
                            continue;

                        FlowType[layer + 1] = 0;
                        flow = ExcessW[layer];
                        K_initial = QK[layer, 0];
                        theta = QTheta[layer, 0];

                        for (int istep = 0; istep <= step; ++istep)
                        {
                            if (FlowType[layer + 1] > 1)
                                break;

                            if (theta > theta_steps[istep, layer])
                            {
                                available = (theta - theta_steps[istep, layer]) * soilPhysical.Thickness[layer];
                                if (available >= K_steps[istep, layer])
                                {
                                    FlowType[layer + 1] = 2;
                                    available = K_steps[istep, layer];
                                }
                                flow += available;
                                theta -= available / soilPhysical.Thickness[layer];
                            }
                        }

                        if (flow > K_initial)
                            flow = K_initial;

                        if (InterfaceFlow[layer + 1] + flow > soilPhysical.KS[layer])
                        {
                            InterfaceFlow[layer + 1] = Math.Max(InterfaceFlow[layer + 1], soilPhysical.KS[layer]);
                            flow = Math.Max(0.0, soilPhysical.KS[layer] - InterfaceFlow[layer + 1]);
                        }
                        else
                        {
                            InterfaceFlow[layer + 1] += flow;
                        }

                        if (flow < ExcessW[layer])
                        {
                            BackFlow[layer] += ExcessW[layer] - flow;
                            if (InterfaceFlow[layer + 1] >= 0.0)
                                FlowType[layer + 1] = 3;
                            else
                                FlowType[layer + 1] = -3;
                        }

                        NewWater[layer] -= flow;
                        ExcessW[layer] = 0;
                        Redistribute(layer, flow, 1);
                        NewWater[layer + 1] += flow;
                        Redistribute(layer + 1, flow, -1);
                        ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                    }

                    if (FlowType[num_layers] == 0)
                    {
                        flow = ExcessW[num_layers - 1];
                        K_initial = QK[num_layers - 1, 0];
                        theta = QTheta[num_layers - 1, 0];

                        for (int istep = 0; istep <= step; ++istep)
                        {
                            if (FlowType[num_layers] > 1)
                                break;

                            if (theta > theta_steps[istep, num_layers - 1])
                            {
                                available = (theta - theta_steps[istep, num_layers - 1]) * soilPhysical.Thickness[num_layers - 1];
                                if (available >= K_steps[istep, num_layers - 1])
                                {
                                    FlowType[num_layers] = 2;
                                    available = K_steps[istep, num_layers - 1];
                                }
                                flow += available;
                                theta -= available / soilPhysical.Thickness[num_layers - 1];
                            }
                        }

                        if (flow > K_initial)
                            flow = K_initial;

                        if (InterfaceFlow[num_layers] + flow > soilPhysical.KS[num_layers - 1])
                        {
                            InterfaceFlow[num_layers] = Math.Max(InterfaceFlow[num_layers], soilPhysical.KS[num_layers - 1]);
                            flow = Math.Max(0.0, soilPhysical.KS[num_layers - 1] - InterfaceFlow[num_layers]);
                        }
                        else
                        {
                            InterfaceFlow[num_layers] += flow;
                        }

                        if (flow < ExcessW[num_layers - 1])
                        {
                            BackFlow[num_layers - 1] += ExcessW[num_layers - 1] - flow;
                        }

                        NewWater[num_layers - 1] -= flow;
                        ExcessW[num_layers - 1] = 0.0;
                        Redistribute(num_layers - 1, flow, 1);
                        FlowType[num_layers] = 0;
                    }
                }
            }

            // TODO: missing backup flow

            #endregion


            #region Gravitational flow (not in use)

            //for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    ExcessW[layer] = 0.0;
            //    BackFlow[layer] = 0.0;
            //}

            //MaxIteration = 10;
            //for (int it = 0; it < MaxIteration; ++it)
            //{
            //    double aWater;
            //    double pFlow;

            //    for (int layer = 0; layer < num_layers - 1; ++layer)
            //    {
            //        if (FlowType[layer + 1] > 1)
            //            continue;

            //        aWater = AvailableWater(layer);
            //        pFlow = aWater + ExcessW[layer];
            //        if (MathUtilities.FloatsAreEqual(0.0, pFlow, 0.1))
            //        {
            //            if (InterfaceFlow[layer + 1] + pFlow >= 0)
            //                FlowType[layer + 1] = 3;
            //            else
            //                System.Console.WriteLine("Flow should not be negative at this stage.");
            //            continue;
            //        }

            //        dist = QDepth[layer + 1, 0] - QDepth[layer, num_Qpoints - 1];
            //        psi_diff = Qpsi[layer, num_Qpoints - 1] - Qpsi[layer + 1, 0];
            //        PFlow[layer + 1] = (QK[layer, num_Qpoints - 1] + QK[layer + 1, 0]) / 2 * (psi_diff / dist + 1.0);
            //        PFlow[layer + 1] = Math.Max(0.0, PFlow[layer + 1]);
            //        MaxFlow = (soilPhysical.KS[layer] + soilPhysical.KS[layer + 1]) / 2;
            //        PFlow[layer + 1] = Math.Min(PFlow[layer + 1], MaxFlow);

            //        aWater = Math.Min(aWater, PFlow[layer + 1]) / 2.0;
            //        pFlow = ExcessW[layer] + aWater;

            //        if (pFlow + InterfaceFlow[layer + 1] > MaxFlow)
            //        {
            //            pFlow = Math.Max(0.0, MaxFlow - InterfaceFlow[layer + 1]);
            //            if (pFlow < ExcessW[layer])
            //            {
            //                BackFlow[layer] += ExcessW[layer] - pFlow;
            //            }

            //            FlowType[layer + 1] = 2;
            //        }

            //        InterfaceFlow[layer + 1] += pFlow;
            //        NewWater[layer] -= pFlow;
            //        NewWater[layer + 1] += pFlow;
            //        if (pFlow > ExcessW[layer])
            //            Redistribute(layer, pFlow - ExcessW[layer], 1);
            //        ExcessW[layer] = Math.Max(0.0, ExcessW[layer] - pFlow);
            //        if (ExcessW[layer + 1] > 0.0)
            //        {
            //            ExcessW[layer + 1] += pFlow;
            //        }
            //        else
            //        {
            //            // is_Saturated is not used in the if () as it wasn't updated after drainage.
            //            is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1] / soilPhysical.SATmm[layer + 1], 0.9999);
            //            ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //        }
            //        Redistribute(layer + 1, pFlow - ExcessW[layer], -1);
            //    }

            //    if (FlowType[num_layers] != 2)
            //    {
            //        aWater = AvailableWater(num_layers - 1);
            //        pFlow = aWater + ExcessW[num_layers - 1];
            //        PFlow[num_layers] = QK[num_layers - 1, num_Qpoints - 1];
            //        MaxFlow = soilPhysical.KS[num_layers - 1];
            //        aWater = Math.Min(aWater, PFlow[num_layers]) / 2.0;
            //        pFlow = aWater + ExcessW[num_layers - 1];
            //        if (pFlow + InterfaceFlow[num_layers] > MaxFlow)
            //        {
            //            pFlow = Math.Max(0.0, MaxFlow - InterfaceFlow[num_layers]);
            //            if (pFlow < ExcessW[num_layers - 1])
            //            {
            //                BackFlow[num_layers - 1] += ExcessW[num_layers - 1] - pFlow;
            //            }
            //            FlowType[num_layers] = 2;
            //        }

            //        InterfaceFlow[num_layers] += pFlow;
            //        NewWater[num_layers - 1] -= pFlow;
            //        if (pFlow > ExcessW[num_layers - 1])
            //            Redistribute(num_layers - 1, pFlow - ExcessW[num_layers - 1], 1);
            //        ExcessW[num_layers - 1] = Math.Max(0.0, ExcessW[num_layers - 1] - pFlow);
            //    }
            //}

            //for (int layer = 1; layer < num_layers; ++layer)
            //{
            //    if (FlowType[layer] == 0)
            //    {
            //        if (InterfaceFlow[layer] > 0.0)
            //            FlowType[layer] = 1;
            //    }
            //    else if (FlowType[layer] == 3)
            //    {
            //        if (InterfaceFlow[layer] > 0.0)
            //            FlowType[layer] = 1;
            //        else
            //            FlowType[layer] = 0;
            //    }
            //}

            //// Backup flow
            //for (int layer = num_layers - 1; layer > 0; --layer)
            //{
            //    if (BackFlow[layer] > 0.0)
            //    {
            //        if (FlowType[layer] > 0)
            //        {
            //            InterfaceFlow[layer] -= BackFlow[layer];

            //            if (InterfaceFlow[layer] < 0.0)
            //                FlowType[layer] = -2;
            //            else
            //                FlowType[layer] = 2;
            //        }
            //        else
            //        {
            //            InterfaceFlow[layer] = -BackFlow[layer];
            //            FlowType[layer] = -1;
            //        }

            //        NewWater[layer] -= BackFlow[layer];
            //        ExcessW[layer] -= BackFlow[layer];
            //        NewWater[layer - 1] += BackFlow[layer];

            //        if (BackFlow[layer - 1] > 0.0)
            //        {
            //            BackFlow[layer - 1] += BackFlow[layer];
            //            ExcessW[layer - 1] += BackFlow[layer];
            //        }
            //        else
            //        {
            //            Redistribute(layer - 1, -BackFlow[layer], 1);
            //            is_Saturated[layer - 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer - 1] / soilPhysical.SATmm[layer - 1], 0.9999);
            //            BackFlow[layer - 1] = Math.Max(0.0, NewWater[layer - 1] - soilPhysical.SATmm[layer - 1]);
            //            ExcessW[layer - 1] = BackFlow[layer - 1];
            //        }
            //    }
            //}

            //if (BackFlow[0] > 0.0)
            //{
            //    InterfaceFlow[0] -= BackFlow[0];

            //    NewWater[0] -= BackFlow[0];
            //    ExcessW[0] -= BackFlow[0];
            //}

            //// Excess water should be 0 at the end.
            //for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    if (!MathUtilities.FloatsAreEqual(ExcessW[layer], 0.0, 1e-6))
            //        System.Console.WriteLine("Excess water is not 0 after forced flow.");
            //}

            #endregion


            #region Matrix flow

            for (int layer = 0; layer < num_layers - 1; ++layer)
            {
                double matrixFlow;

                if (Math.Abs(FlowType[layer + 1]) > 2)
                    continue;

                matrixFlow = UnsatFlow(layer);
                MaxFlow = (soilPhysical.KS[layer] + soilPhysical.KS[layer + 1]) / 2.0;
                if (matrixFlow > 0.0)
                {
                    switch (FlowType[layer + 1])
                    {
                        case -1:
                            continue;

                        case 0:
                            FlowType[layer + 1] = 1;
                            break;

                        default:
                            break;
                    }

                    if (matrixFlow + InterfaceFlow[layer + 1] >= MaxFlow)
                    {
                        FlowType[layer + 1] = 2;
                        matrixFlow = MaxFlow - InterfaceFlow[layer + 1];
                        matrixFlow = Math.Max(0.0, matrixFlow);
                    }
                }
                else if (matrixFlow < 0.0)
                {
                    switch (FlowType[layer + 1])
                    {
                        case -1:
                            break;

                        case 0:
                            FlowType[layer + 1] = -1;
                            break;

                        case 1:
                            continue;
                    }

                    if (matrixFlow + InterfaceFlow[layer + 1] <= -MaxFlow)
                    {
                        FlowType[layer + 1] = -2;
                        matrixFlow = MaxFlow - InterfaceFlow[layer + 1];
                        matrixFlow = Math.Min(0.0, matrixFlow);
                    }
                }

                InterfaceFlow[layer + 1] += matrixFlow;
                NewWater[layer] -= matrixFlow;
                Redistribute(layer, matrixFlow, 1);
                NewWater[layer + 1] += matrixFlow;
                Redistribute(layer + 1, matrixFlow, -1);
            }

            #endregion

            #region Gravitational and matrix flow?

            //// Calculate flow using the initial status, which will be the max flow for that interface (within and between layers).
            //for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    // Calculate initial flow within a layer.
            //    Redistribute(layer, -1);

            //    ExcessW[layer] = 0.0;
            //    BackFlow[layer] = 0.0;
            //}

            //for (int layer = 0; layer < num_layers - 1; ++layer)
            //{
            //    dist = QDepth[layer + 1, 0] - QDepth[layer, num_Qpoints - 1];
            //    psi_diff = Qpsi[layer, num_Qpoints - 1] - Qpsi[layer + 1, 0];
            //    PFlow[layer + 1] = (QK[layer, num_Qpoints - 1] + QK[layer + 1, 0]) / 2 * (psi_diff / dist + 1.0);

            //    //if (!is_Saturated[layer + 1])
            //    //{
            //    //    dist = QDepth[layer + 1, 0] - QDepth[layer, num_Qpoints - 1];
            //    //    psi_diff = Qpsi[layer, num_Qpoints - 1] - Qpsi[layer + 1, 0];
            //    //    PFlow[layer] = (QK[layer, num_Qpoints - 1] + QK[layer + 1, 0]) / 2 * (psi_diff / dist + 1.0);
            //    //}
            //    //else
            //    //{
            //    //    dist = soilPhysical.ThicknessCumulative[layer] - QDepth[layer, num_Qpoints - 1];
            //    //    psi_diff = Qpsi[layer, num_Qpoints - 1] - 0.0;
            //    //    PFlow[layer] = QK[layer, num_Qpoints - 1] * (psi_diff / dist + 1.0);
            //    //}

            //    if (PFlow[layer + 1] < 0.0)
            //    {
            //        MaxFlow = -(QK[layer, num_Qpoints - 1] + QK[layer + 1, 0]) / 2;
            //        PFlow[layer + 1] = Math.Max(PFlow[layer + 1], MaxFlow);
            //        if (FlowType[layer + 1] == 0)
            //            FlowType[layer + 1] = -1;
            //    }
            //    else
            //    {
            //        MaxFlow = QK[layer, num_Qpoints - 1] + QK[layer + 1, 0];
            //        PFlow[layer + 1] = Math.Min(PFlow[layer + 1], MaxFlow);
            //        if (FlowType[layer + 1] == 0)
            //            FlowType[layer + 1] = 1;
            //    }
            //}

            //// TODO: other bottom boundary types.
            //PFlow[num_layers] = QK[num_layers - 1, num_Qpoints - 1];


            //// Iteration
            //// Move water until PFlow is reached or other criteria are met.
            //int i = 0;
            //for (; i < MaxIteration; ++i)
            //{
            //    bool exitCon = true;
            //    for (int layer = 0; layer < num_layers - 1; ++layer)
            //    {
            //        switch (FlowType[layer + 1])
            //        {
            //            case -2:
            //                break;

            //            case -1:
            //                UpFlow = UnsatFlow(layer);
            //                if (ExcessW[layer] > 0.0)
            //                {
            //                    FlowType[layer + 1] = -2;
            //                    BackFlow[layer] += ExcessW[layer];
            //                    break;
            //                }
            //                if (UpFlow > 0.0)
            //                {
            //                    FlowType[layer + 1] = -2;
            //                }
            //                else
            //                {
            //                    if (InterfaceFlow[layer + 1] + UpFlow > PFlow[layer + 1])
            //                    {
            //                        InterfaceFlow[layer + 1] += UpFlow;
            //                        if (Math.Abs(UpFlow) < 1e-6)
            //                        {
            //                            FlowType[layer + 1] = -2;
            //                        }
            //                        else
            //                        {
            //                            exitCon = false;
            //                        }
            //                    }
            //                    else
            //                    {
            //                        FlowType[layer + 1] = -2;
            //                        UpFlow = PFlow[layer + 1] - InterfaceFlow[layer + 1];
            //                        InterfaceFlow[layer + 1] = PFlow[layer + 1];
            //                    }

            //                    NewWater[layer] -= UpFlow;
            //                    Redistribute(layer, UpFlow, 1);
            //                    NewWater[layer + 1] += UpFlow;
            //                    Redistribute(layer + 1, UpFlow, -1);
            //                    // NewWater[layer] shouldn't reach saturation here.
            //                }
                            
            //                break;

            //            case 0:
            //                System.Console.WriteLine("FlowType shouldn't be 0 here.");
            //                break;

            //            case 1:
            //                DownFlow = UnsatFlow(layer);
            //                if (DownFlow < 0.0)
            //                {
            //                    FlowType[layer + 1] = 2;
            //                    BackFlow[layer] += ExcessW[layer];
            //                    // ExcessW shouldn't be greater than 0 here.
            //                }
            //                else
            //                {
            //                    DownFlow += ExcessW[layer];
            //                    if (InterfaceFlow[layer + 1] + DownFlow < PFlow[layer + 1])
            //                    {
            //                        InterfaceFlow[layer + 1] += DownFlow;
            //                        Redistribute(layer, DownFlow - ExcessW[layer], 1);
            //                        exitCon = false;
            //                    }
            //                    else
            //                    {
            //                        FlowType[layer + 1] = 2;
            //                        DownFlow = PFlow[layer + 1] - InterfaceFlow[layer + 1];
            //                        InterfaceFlow[layer + 1] = PFlow[layer + 1];
            //                        if (DownFlow < ExcessW[layer])
            //                            BackFlow[layer] += ExcessW[layer] - DownFlow;
            //                        else
            //                            Redistribute(layer, DownFlow - ExcessW[layer], 1);
            //                    }

            //                    NewWater[layer] -= DownFlow;
            //                    NewWater[layer + 1] += DownFlow;
            //                    Redistribute(layer + 1, DownFlow, -1);
            //                    ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //                }

            //                break;
            //            case 2:
            //                break;
            //            default:
            //                System.Console.WriteLine("Error: unkown flow type.");
            //                break;
            //        }
            //    }

            //    // Last layer
            //    switch (FlowType[num_layers - 1])
            //    {
            //        case 0:
            //        case 1:
            //            DownFlow = UnsatFlow(num_layers - 1);
            //            DownFlow += ExcessW[num_layers - 1];
            //            if (MathUtilities.FloatsAreEqual(DownFlow, 0, 1e-6))
            //            {
            //                FlowType[num_layers] = 2;
            //                InterfaceFlow[num_layers] += DownFlow;
            //            }
            //            else
            //            {
            //                if (InterfaceFlow[num_layers] + DownFlow < PFlow[num_layers])
            //                {
            //                    InterfaceFlow[num_layers] += DownFlow;
            //                    exitCon = false;
            //                }
            //                else
            //                {
            //                    FlowType[num_layers] = 2;
            //                    DownFlow = PFlow[num_layers] - InterfaceFlow[num_layers];
            //                    InterfaceFlow[num_layers] = PFlow[num_layers];
            //                    if (DownFlow < ExcessW[num_layers - 1])
            //                        BackFlow[num_layers - 1] += ExcessW[num_layers - 1] - DownFlow;
            //                }
            //            }
            //            NewWater[num_layers - 1] -= DownFlow;
            //            Redistribute(num_layers - 1, DownFlow - ExcessW[num_layers - 1], 1);

            //            break;
            //        case 2:
            //            break;

            //    }

            //    // Backup flow
            //    for (int layer = num_layers - 1; layer > 0; --layer)
            //    {
            //        if (BackFlow[layer] > 0.0)
            //        {
            //            InterfaceFlow[layer] -= BackFlow[layer];

            //            if (FlowType[layer] > 0)
            //            {
            //                if (InterfaceFlow[layer] < 0.0)
            //                    FlowType[layer] = -2;
            //                else
            //                    FlowType[layer] = 2;
            //            }

            //            NewWater[layer] -= BackFlow[layer];
            //            ExcessW[layer] -= BackFlow[layer];
            //            NewWater[layer - 1] += BackFlow[layer];

            //            if (BackFlow[layer - 1] > 0.0)
            //            {
            //                BackFlow[layer - 1] += BackFlow[layer];
            //                ExcessW[layer - 1] += BackFlow[layer];
            //            }
            //            else
            //            {
            //                Redistribute(layer - 1, -BackFlow[layer], 1);
            //                is_Saturated[layer - 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer - 1] / soilPhysical.SATmm[layer - 1], 0.9999);
            //                BackFlow[layer - 1] = Math.Max(0.0, NewWater[layer - 1] - soilPhysical.SATmm[layer - 1]);
            //                ExcessW[layer - 1] = BackFlow[layer - 1];
            //            }
            //        }
            //    }

            //    if (BackFlow[0] > 0.0)
            //    {
            //        InterfaceFlow[0] -= BackFlow[0];

            //        NewWater[0] -= BackFlow[0];
            //        ExcessW[0] -= BackFlow[0];
            //    }

            //    // Water flow within each layer.
            //    for (int layer = 0; layer < num_layers; ++layer)
            //    {
            //        Redistribute(layer, 1);
            //    }

            //    if (exitCon == true)
            //        break;
            //}
            //System.Console.WriteLine("Iteration ended when i euqals to: " + i + ".");
            #endregion

            #region Not in use
            //// add evaporation
            //if (Infiltration > 0.0)
            //{
            //    // TODO: This flux might limit infiltration when soil is relatively dry.
            //    DownFlow = QK[0, 0] * ((Pond - (Qpsi[0, 0]) / QDepth[0, 0] + 1));
            //    MaxFlow = Math.Max(DownFlow, 2 * soilPhysical.KS[0]);
            //    InterfaceFlow[0] = Math.Min(MaxFlow, Infiltration);
            //    NewWater[0] += InterfaceFlow[0];
            //    Redistribute(0, InterfaceFlow[0], -1);
            //    if (ExcessW[0] > 0.0)
            //        ExcessW[0] += InterfaceFlow[0];
            //    else
            //    {
            //        is_Saturated[0] = MathUtilities.IsGreaterThanOrEqual(NewWater[0], soilPhysical.SATmm[0]);
            //        ExcessW[0] = Math.Max(0.0, NewWater[0] - soilPhysical.SATmm[0]);
            //    }
            //}

            //for (int layer = 0; layer < num_layers - 1; ++layer)
            //{
            //    if (FlowType[layer] == 1)
            //    {
            //        if (ExcessW[layer] > 0.0)
            //        {
            //            FlowType[layer + 1] = 1;
            //            InterfaceFlow[layer + 1] = ExcessW[layer];
            //            MaxFlow = 2 * soilPhysical.KS[layer + 1];
            //            BackFlow[layer] = InterfaceFlow[layer + 1] - MaxFlow;
            //            if (BackFlow[layer] > 0.0)
            //            {
            //                InterfaceFlow[layer + 1] -= BackFlow[layer];
            //                Redistribute(layer + 1, InterfaceFlow[layer + 1], -1);
            //            }
            //            else
            //            {
            //                Redistribute(layer + 1, InterfaceFlow[layer + 1], -1);
            //                DownFlow = UnsatFlow(layer);
            //                if (MathUtilities.IsGreaterThanOrEqual(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
            //                    DownFlow = MaxFlow - InterfaceFlow[layer + 1];

            //                InterfaceFlow[layer + 1] += DownFlow;
            //                Redistribute(layer, DownFlow, 1);
            //                Redistribute(layer + 1, DownFlow, -1);
            //            }
            //            NewWater[layer] -= InterfaceFlow[layer + 1];
            //            NewWater[layer + 1] += InterfaceFlow[layer + 1];
            //            if (ExcessW[layer + 1] > 0.0)
            //                ExcessW[layer + 1] += InterfaceFlow[layer + 1];
            //            else
            //            {
            //                is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
            //                ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //            }
            //        }
            //        else
            //        {
            //            DownFlow = UnsatFlow(layer);
            //            InterfaceFlow[layer + 1] = DownFlow;
            //            if (DownFlow > 0.0)
            //            {
            //                FlowType[layer + 1] = 1;
            //                Redistribute(layer, DownFlow, 1);
            //                NewWater[layer] -= InterfaceFlow[layer + 1];
            //                Redistribute(layer + 1, DownFlow, -1);
            //                NewWater[layer + 1] += InterfaceFlow[layer + 1];
            //                if (ExcessW[layer + 1] > 0.0)
            //                    ExcessW[layer + 1] += InterfaceFlow[layer + 1];
            //                else
            //                {
            //                    is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
            //                    ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //                }
            //            }
            //        }
            //    }
            //    else
            //    {
            //        if (ExcessW[layer] > 0.0)
            //        {
            //            double amount;
            //            bool is_upflow = true;
            //            FlowType[layer + 1] = 1;
            //            InterfaceFlow[layer] = 0.0;
            //            InterfaceFlow[layer + 1] = 0.0;

            //            do
            //            {
            //                if (ExcessW[layer] < 2.0)
            //                    amount = ExcessW[layer];
            //                else
            //                    amount = 2.0;
            //                ExcessW[layer] -= 2.0;

            //                if (is_upflow)
            //                {
            //                    // This layer is fully saturated. The downward flow is calculate assuming pressure building in the whole layer.
            //                    double up = -QK[layer - 1, num_Qpoints - 1] * ((Qpsi[layer - 1, num_Qpoints - 1] - 0.0) / (soilPhysical.ThicknessCumulative[layer - 1] - QDepth[layer - 1, num_Qpoints - 1]) + 1);
            //                    double down = QK[layer + 1, 0] * ((0.0 - Qpsi[layer + 1, 0]) / (QDepth[layer + 1, 0] - soilPhysical.ThicknessCumulative[layer]) + 2);
            //                    if (up > 0)
            //                    {
            //                        UpFlow = up / (up + down) * amount;
            //                        Theta = QTheta[layer - 1, num_Qpoints - 1] + UpFlow / (soilPhysical.Thickness[layer - 1] * QWeight[layer - 1, num_Qpoints - 1]);
            //                        if (MathUtilities.IsGreaterThanOrEqual(Theta, soilPhysical.SAT[layer - 1]))
            //                        {
            //                            UpFlow -= (Theta - soilPhysical.SAT[layer - 1]) * (soilPhysical.Thickness[layer - 1] * QWeight[layer - 1, num_Qpoints - 1]);
            //                            is_upflow = false;
            //                        }
            //                        InterfaceFlow[layer] += -UpFlow;
            //                        Redistribute(layer - 1, -UpFlow, 1);
            //                        //Qpsi[layer - 1, 2] = Suction(layer - 1, QTheta[layer - 1, 2]);
            //                        //QK[layer - 1, 2] = SimpleK(layer - 1, Qpsi[layer - 1, 2]);
            //                        Qpsi[layer - 1, num_Qpoints - 1] = hydraulicModel.get_h(layer - 1, QTheta[layer - 1, num_Qpoints - 1]);
            //                        QK[layer - 1, num_Qpoints - 1] = hydraulicModel.get_K(layer - 1, Qpsi[layer - 1, num_Qpoints - 1]);

            //                    }
            //                    else
            //                    {
            //                        is_upflow = false;
            //                        UpFlow = 0.0;
            //                    }
            //                }
            //                else
            //                    UpFlow = 0.0;

            //                DownFlow = amount - UpFlow;

            //                InterfaceFlow[layer + 1] += DownFlow;
            //                Redistribute(layer + 1, DownFlow, -1);
            //            } while (ExcessW[layer] > 0.0);

            //            // continue unsatflow
            //            if (is_upflow)
            //            {
            //                UpFlow = -UnsatFlow(layer - 1);
            //                if (UpFlow > 0.0)
            //                {
            //                    InterfaceFlow[layer] += -UpFlow;
            //                    Redistribute(layer - 1, -UpFlow, 1);
            //                    Redistribute(layer, -UpFlow, -1);
            //                }
            //            }

            //            MaxFlow = 2 * soilPhysical.KS[layer + 1];
            //            if (MathUtilities.IsLessThan(InterfaceFlow[layer + 1], MaxFlow))
            //            {
            //                DownFlow = UnsatFlow(layer);
            //                if (MathUtilities.IsGreaterThan(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
            //                    DownFlow = MaxFlow - InterfaceFlow[layer + 1];
            //                InterfaceFlow[layer + 1] += DownFlow;
            //                Redistribute(layer, DownFlow, 1);
            //            }

            //            NewWater[layer - 1] -= InterfaceFlow[layer];
            //            NewWater[layer] += InterfaceFlow[layer];
            //            NewWater[layer] -= InterfaceFlow[layer + 1];
            //            NewWater[layer + 1] += InterfaceFlow[layer + 1];
            //            if (ExcessW[layer + 1] > 0.0)
            //                ExcessW[layer + 1] += InterfaceFlow[layer + 1];
            //            else
            //            {
            //                is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
            //                ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //            }
            //        }
            //        else
            //        {
            //            UpFlow = -UnsatFlow(layer - 1);
            //            InterfaceFlow[layer] = -UpFlow;
            //            Redistribute(layer - 1, -UpFlow, 1);
            //            Redistribute(layer, -UpFlow, -1);
            //            NewWater[layer - 1] -= InterfaceFlow[layer];
            //            NewWater[layer] += InterfaceFlow[layer];

            //            DownFlow = UnsatFlow(layer);
            //            InterfaceFlow[layer + 1] = DownFlow;
            //            if (DownFlow > 0.0)
            //            {
            //                FlowType[layer + 1] = 1;
            //                Redistribute(layer, DownFlow, 1);
            //                Redistribute(layer + 1, DownFlow, -1);
            //                NewWater[layer] -= InterfaceFlow[layer + 1];
            //                NewWater[layer + 1] += InterfaceFlow[layer + 1];
            //                if (ExcessW[layer + 1] > 0.0)
            //                    ExcessW[layer + 1] += InterfaceFlow[layer + 1];
            //                else
            //                {
            //                    is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
            //                    ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
            //                }
            //            }
            //        }
            //    }
            //}

            //// Last layer.
            //// TODO: other boundary type.
            //{
            //    int layer = num_layers- 1;
            //    InterfaceFlow[layer + 1] = 0.0;

            //    if (FlowType[layer] == 1)
            //    {
            //        if (ExcessW[layer] > 0.0)
            //        {
            //            InterfaceFlow[layer + 1] = ExcessW[layer];
            //            MaxFlow = soilPhysical.KS[layer];
            //            BackFlow[layer] = InterfaceFlow[layer + 1] - MaxFlow;
            //            if (BackFlow[layer] > 0.0)
            //            {
            //                InterfaceFlow[layer + 1] -= BackFlow[layer];
            //            }
            //            else
            //            {
            //                DownFlow = UnsatFlow(layer);
            //                if (MathUtilities.IsGreaterThan(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
            //                    DownFlow = MaxFlow - InterfaceFlow[layer + 1];
            //                InterfaceFlow[layer + 1] += DownFlow;
            //                Redistribute(layer, DownFlow, 1);
            //            }
            //        }
            //        else
            //        {
            //            DownFlow = UnsatFlow(layer);
            //            InterfaceFlow[layer + 1] = DownFlow;
            //            Redistribute(layer, DownFlow, 1);
            //        }
            //    }
            //    else
            //    {
            //        UpFlow = -UnsatFlow(layer - 1);
            //        InterfaceFlow[layer] = -UpFlow;
            //        Redistribute(layer - 1, -UpFlow, 1);
            //        Redistribute(layer, -UpFlow, -1);
            //        NewWater[layer - 1] -= InterfaceFlow[layer];
            //        NewWater[layer] += InterfaceFlow[layer];

            //        DownFlow = UnsatFlow(layer);
            //        InterfaceFlow[layer + 1] = DownFlow;
            //        Redistribute(layer, DownFlow, 1);
            //    }

            //    NewWater[layer] -= InterfaceFlow[layer + 1];
            //}

            //// Backup flow
            //for (int layer = num_layers - 1; layer > 0; --layer)
            //{
            //    if (MathUtilities.IsGreaterThan(NewWater[layer], soilPhysical.SATmm[layer]))
            //    {
            //        BackFlow[layer] = NewWater[layer] - soilPhysical.SATmm[layer];
            //        InterfaceFlow[layer] -= BackFlow[layer];
            //        Redistribute(layer - 1, -BackFlow[layer], 1);
            //        NewWater[layer - 1] += BackFlow[layer];
            //        NewWater[layer] -= BackFlow[layer];
            //    }
            //}
            //if (MathUtilities.IsGreaterThan(NewWater[0], soilPhysical.SATmm[0]))
            //{
            //    BackFlow[0] = NewWater[0] - soilPhysical.SATmm[0];
            //    InterfaceFlow[0] -= BackFlow[0];
            //}

            //// Unsaturated flow within each layer.
            //for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    Redistribute(layer, 1);
            //}

            #endregion
        }

        /// <summary>
        /// Redistribute water within one layer.
        /// </summary>
        /// <param name="layer"></param>
        /// <param name="method"></param>
        private void Redistribute(int layer, int method = 0)
        {
            // This method is used to calculate a)water distribution after changes outside of this class (root uptake),
            // and b) redistribution at the end.

            const double eff = 1;
            double[] depth = new double[num_Qpoints];

            // Mass balance
            double TotalWater_old = 0;
            double TotalWater_new = 0;

            switch (method)
            {
                case 0:
                    // To adjust water content after root uptake.

                    double uptake;
                    double[] PAW = new double[num_Qpoints];

                    if (!MathUtilities.FloatsAreEqual(OldWater[layer], Water[layer]))
                    {
                        // Subtract water evenly from each section unless not enough water is available in some sections.
                        // TODO: some sections may not have enough plant available water.
                        uptake = OldWater[layer] - Water[layer];

                        for (int n = 0; n <= num_Qpoints; ++n)
                        {
                            depth[n] = soilPhysical.Thickness[layer] * QWeight[layer, n];
                            QTheta[layer, n] -= uptake * QWeight[layer, n] / depth[n];
                        }

                        for (int n = 0; n <= num_Qpoints; ++n)
                        {
                            Qpsi[layer, n] = hydraulicModel.get_h(layer, QTheta[layer, n]);
                            QK[layer, n] = hydraulicModel.get_K(layer, Qpsi[layer, n]);
                        }
                    }

                    break;

                case -1:
                    // Calculate max flux between sublayers.
                    // These values are set for case 1.
                    //for (int n = 0; n < num_Qpoints - 1; ++n)
                    //{
                    //    QFlux[layer, n] = (QK[layer, n] + QK[layer, n + 1]) / 2
                    //            * ((Qpsi[layer, n] - Qpsi[layer, n + 1]) / (QDepth[layer, n + 1] - QDepth[layer, n]) + 1);
                    //    QFlux[layer, n] = Math.Min(QFlux[layer, n], QK[layer, n] + QK[layer, n + 1]);
                    //}

                    break;

                case 1:
                    // TODO: not modified for new method. Not in use.
                    // Internal unsaturated flow after calculating all interface flow.
                    // Scan the profile to get the distribution type first.
                    // Only three points case considered right now.
                    double theta_mean = 0;
                    double psi_mean;
                    double cap;
                    double delta;

                    double[] theta_eq = new double[num_Qpoints];
                    double[] water_diff = new double[num_Qpoints];

                    double[] flow = new double[num_Qpoints - 1];

                    for (int n = 0; n < num_Qpoints; ++n)
                    {
                        depth[n] = soilPhysical.Thickness[layer] * QWeight[layer, n];
                        theta_mean += QTheta[layer, n] * QWeight[layer, n];
                        TotalWater_old += QTheta[layer, n] * depth[n];
                    }

                    // TODO: a new method to accelerate water movement when the soil is very dry.
                    // Hint: the gradient would be much great if the depth is reduced.
                    if (num_Qpoints == 3)
                    {
                        for (int n = 0; n < num_Qpoints - 1; ++n)
                        {
                            flow[n] = (QK[layer, n] + QK[layer, n + 1]) / 2
                                * ((Qpsi[layer, n] - Qpsi[layer, n + 1]) / (QDepth[layer, n + 1] - QDepth[layer, n]) + 1);

                            // flow[n] /= 2.0;
                        }

                        if (MathUtilities.FloatsAreEqual(theta_mean, soilPhysical.SAT[layer]))
                            break;

                        // Calculate the equilibrium distribution.
                        delta = (depth[0] + depth[1]) / 2.0;
                        //psi_mean = Suction(layer, theta_mean);
                        //cap = (SimpleTheta(layer, psi_mean + delta) - SimpleTheta(layer, psi_mean - delta)) / (2 * delta);
                        psi_mean = hydraulicModel.get_h(layer, theta_mean);
                        cap = (hydraulicModel.get_theta(layer, psi_mean + delta) - hydraulicModel.get_theta(layer, psi_mean - delta)) / (2 * delta);
                        theta_eq[0] = theta_mean - cap * delta;
                        theta_eq[1] = theta_mean;
                        theta_eq[2] = theta_mean + cap * delta;
                        if (theta_eq[2] > soilPhysical.SAT[layer])
                        {
                            theta_eq[1] += (theta_eq[2] - soilPhysical.SAT[layer]) * depth[2] / depth[1];
                            theta_eq[2] = soilPhysical.SAT[layer];
                            if (theta_eq[1] > soilPhysical.SAT[layer])
                            {
                                theta_eq[0] += (theta_eq[1] - soilPhysical.SAT[layer]) * depth[1] / depth[0];
                                theta_eq[1] = soilPhysical.SAT[layer];
                            }
                        }
                        for (int n = 0; n < num_Qpoints; ++n)
                        {
                            water_diff[n] = (QTheta[layer, n] - theta_eq[n]) * depth[n];
                        }

                        // TODO: need a better method to determine how fast it is reaching equilibrium.
                        if (water_diff[0] * flow[0] > 0.0)
                        {
                            if (Math.Abs(flow[0]) < Math.Abs(water_diff[0]))
                                QTheta[layer, 0] -= eff * flow[0] / depth[0];
                            else
                                QTheta[layer, 0] -= eff * water_diff[0] / depth[0];
                        }

                        if (water_diff[2] * flow[1] < 0.0)
                        {
                            if (Math.Abs(flow[1]) < Math.Abs(water_diff[2]))
                                QTheta[layer, 2] += eff * flow[1] / depth[2];
                            else
                                QTheta[layer, 2] -= eff * water_diff[2] / depth[2];
                        }

                        water_diff[0] = (QTheta[layer, 0] - theta_eq[0]) * depth[0];
                        water_diff[2] = (QTheta[layer, 2] - theta_eq[2]) * depth[2];
                        water_diff[1] = -(water_diff[0] + water_diff[2]);
                        QTheta[layer, 1] = theta_eq[1] + water_diff[1] / depth[1];
                    }
                    

                    for (int n = 0; n < num_Qpoints; ++n)
                    {
                        //Qpsi[layer, n] = Suction(layer, QTheta[layer, n]);
                        //QK[layer, n] = SimpleK(layer, Qpsi[layer, n]);
                        Qpsi[layer, n] = hydraulicModel.get_h(layer, QTheta[layer, n]);
                        QK[layer, n] = hydraulicModel.get_K(layer, Qpsi[layer, n]);
                        TotalWater_new += QTheta[layer, n] * depth[n];
                    }

                    if (!MathUtilities.FloatsAreEqual(TotalWater_new, TotalWater_old, 1e-6))
                        System.Console.WriteLine("Water mass balance error.");

                    break;

                default:
                    break;
            }
        }

        /// <summary>
        /// Redistribute water flow in and out of one layer.
        /// Inflow is considered as "forced flow" while outflow will only drain the corresponding layer.
        /// </summary>
        /// <param name="layer"></param>
        /// <param name="flow"></param>
        /// <param name="location"></param>
        private void Redistribute(int layer, double flow, int location = 1)
        {
            // Use location to indicate if flow is from top or bottom.
            // 1: water move in or out from the bottom;
            // 0: water move in or out from the middle;
            // -1: water move in or out from the top.
            // Positive flow is inflow for 0 and -1, outflow for 1.

            double available;
            double deficit;
            double meanTheta;
            double meanPsi;
            double meanK;

            meanTheta = NewWater[layer] / soilPhysical.Thickness[layer];
            meanTheta = Math.Min(meanTheta, soilPhysical.SAT[layer]);
            meanPsi = hydraulicModel.get_h(layer, meanTheta);
            meanK = hydraulicModel.get_K(layer, meanPsi);

            QTheta[layer, 0] = meanTheta;

            switch (location)
            {
                case 1:
                    if (flow >= 0.0)
                    {
                        // Flow out from the bottom of this layer.
                        // Use the available DUL water first (from top), then the water above LL15 in last section of the layer.

                        // TODO: using average value for all sections temporarily.
                        for (int q = 0; q <= num_Qpoints; ++q)
                        {
                            QTheta[layer, q] = meanTheta;
                            Qpsi[layer, q] = meanPsi;
                            QK[layer, q] = meanK;
                        }

                        //available = AvailableWater(layer);
                        //if (flow <= available)
                        //{
                        //    for (int q = 0; q < num_Qpoints; ++q)
                        //    {
                        //        available = (QTheta[layer, q] - soilPhysical.DUL[layer]) * (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                        //        if (flow <= available)
                        //        {
                        //            QTheta[layer, q] -= flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                        //            flow = 0;
                        //            break;
                        //        }
                        //        else
                        //        {
                        //            QTheta[layer, q] = soilPhysical.DUL[layer];
                        //            flow -= available;
                        //        }
                        //    }
                        //}
                        //else
                        //{
                        //    for (int q = 0; q < num_Qpoints; ++q)
                        //    {
                        //        QTheta[layer, q] = soilPhysical.DUL[layer];
                        //    }

                        //    flow -= available;

                        //    QTheta[layer, num_Qpoints - 1] -= flow / (QWeight[layer, num_Qpoints - 1] * soilPhysical.Thickness[layer]);
                        //    if (MathUtilities.IsLessThan(QTheta[layer, num_Qpoints - 1], soilPhysical.LL15[layer]))
                        //        System.Console.WriteLine("Too much water flow out from this section.");
                        //}
                    }
                    else
                    {
                        // Water enter the layer from the bottom.
                        // Saturated the lower section before move water up.

                        QTheta[layer, 0] = meanTheta;
                        for (int q = 1; q <= num_Qpoints; ++q)
                        {
                            deficit = (soilPhysical.SAT[layer] - QTheta[layer, q]) * (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                            if (-flow <= deficit)
                            {
                                QTheta[layer, q] -= flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                flow = 0.0;
                                break;
                            }
                            else
                            {
                                QTheta[layer, q] = soilPhysical.SAT[layer];
                                flow += deficit;
                            }
                        }
                    }

                    for (int q = 0; q <= num_Qpoints; ++q)
                    {
                        Qpsi[layer, q] = hydraulicModel.get_h(layer, QTheta[layer, q]);
                        QK[layer, q] = hydraulicModel.get_K(layer, Qpsi[layer, q]);
                    }

                    break;

                case 0:
                    // Evenly distribute in the profile if DUL is not reached.
                    // Extra water above DUL goes to the lower section first.

                    deficit = 0.0;
                    if (flow > 0.0)
                    {
                        for (int q = 1; q <= num_Qpoints; ++q)
                        {
                            deficit += (soilPhysical.DUL[layer] - QTheta[layer, q]) * (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                        }

                        if (flow <= deficit)
                        {
                            for (int q = 1; q <= num_Qpoints; ++q)
                            {
                                deficit = (soilPhysical.DUL[layer] - QTheta[layer, q]) * (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                if (flow <= deficit)
                                {
                                    QTheta[layer, q] += flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                    flow = 0.0;
                                    break;
                                }
                                else
                                {
                                    QTheta[layer, q] = soilPhysical.DUL[layer];
                                    flow -= deficit;
                                }
                            }
                        }
                        else
                        {
                            for (int q = 1; q <= num_Qpoints; ++q)
                            {
                                QTheta[layer, q] = soilPhysical.DUL[layer];
                            }

                            flow -= deficit;
                        }

                        if (flow > 0.0)
                        {
                            for (int q = num_Qpoints - 1; q > 0; --q)
                            {
                                deficit = (soilPhysical.SAT[layer] - QTheta[layer, q]) * (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                if (flow <= deficit)
                                {
                                    QTheta[layer, q] += flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                    flow = 0.0;
                                    break;
                                }
                                else
                                {
                                    QTheta[layer, q] = soilPhysical.DUL[layer];
                                    flow -= deficit;
                                }
                            }
                        }
                    }
                    else
                    {
                        System.Console.WriteLine("Use Distribute without location if water is extracted from the center.");
                    }

                    for (int q = 0; q <= num_Qpoints; ++q)
                    {
                        Qpsi[layer, q] = hydraulicModel.get_h(layer, QTheta[layer, q]);
                        QK[layer, q] = hydraulicModel.get_K(layer, Qpsi[layer, q]);
                    }

                    break;

                case -1:
                    Qpsi[layer, 0] = hydraulicModel.get_h(layer, QTheta[layer, 0]);
                    QK[layer, 0] = hydraulicModel.get_K(layer, Qpsi[layer, 0]);
                    if (flow >= 0.0)
                    {
                        // Fill the layer to DUL from the top, then fill up the layer from bottom section.
                        
                        for (int q = 1; q <= num_Qpoints; ++q)
                        {
                            if (QTheta[layer, q] < soilPhysical.DUL[layer])
                            {
                                deficit = (soilPhysical.DUL[layer] - QTheta[layer, q]) * QWeight[layer, q] * soilPhysical.Thickness[layer];
                                if (flow <= deficit)
                                {
                                    QTheta[layer, q] += flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                    flow = 0.0;
                                    break;
                                }
                                else
                                {
                                    QTheta[layer, q] = soilPhysical.DUL[layer];
                                    flow -= deficit;
                                }
                            }
                        }

                        if (flow > 0.0)
                        {
                            for (int q = 1; q <= num_Qpoints; ++q)
                            {
                                QTheta[layer, q] = meanTheta;
                                Qpsi[layer, q] = meanPsi;
                                QK[layer, q] = meanK;
                            }

                            for (int q = num_Qpoints - 1; q > 0; --q)
                            {
                                //deficit = (soilPhysical.SAT[layer] - QTheta[layer, q]) * QWeight[layer, q] * soilPhysical.Thickness[layer];
                                //if (flow > deficit)
                                //{
                                //    QTheta[layer, q] = soilPhysical.SAT[layer];
                                //    flow -= deficit;
                                //}
                                //else
                                //{
                                //    QTheta[layer, q] += flow / (QWeight[layer, q] * soilPhysical.Thickness[layer]);
                                //    flow = 0.0;
                                //}
                            }
                        }

                        for (int q = 0; q <= num_Qpoints; ++q)
                        {
                            Qpsi[layer, q] = hydraulicModel.get_h(layer, QTheta[layer, q]);
                            QK[layer, q] = hydraulicModel.get_K(layer, Qpsi[layer, q]);
                        }

                    }
                    else
                    {
                        // Only the water in up section is available for upward flow.
                        
                        available = (QTheta[layer, 1] - soilPhysical.LL15[layer]) * QWeight[layer, 1] * soilPhysical.Thickness[layer];
                        available = Math.Max(0.0, available);
                        if (flow <= available)
                        {
                            QTheta[layer, 1] -= flow / (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                            Qpsi[layer, 1] = hydraulicModel.get_h(layer, QTheta[layer, 1]);
                            QK[layer, 1] = hydraulicModel.get_K(layer, Qpsi[layer, 1]);
                        }
                        else
                        {
                            System.Console.WriteLine("Too much water flow out from the top of this layer.");
                        }
                    }

                    break;

                default:
                    break;
            }
        }

        /// <summary>
        /// Calculate the unsaturated flow between this layer and the layer below.
        /// </summary>
        /// <param name="layer"></param>
        /// <returns></returns>
        private double UnsatFlow(int layer)
        {
            // Calculate unsaturated flow between this layer and the one below.

            const double tolerance = 1e-9;
            double water_a;
            double water_d;
            double flow;
            double old_flow;
            double max_flow;
            double min_flow;
            double grad;
            double psi_diff;
            double depth;

            double theta_1;
            double theta_2;
            double psi_1;
            double psi_2;

            if (layer != num_layers - 1)
            {
                depth = (QDepth[layer + 1, 1] - QDepth[layer, num_Qpoints]);
                psi_diff = (Qpsi[layer, num_Qpoints] - Qpsi[layer + 1, 1]);
                grad = psi_diff / depth + 1.0;
                flow = (QK[layer, num_Qpoints] + QK[layer + 1, 1]) / 2 * grad;
                if (Math.Abs(flow) < tolerance)
                    return flow;
                min_flow = 0.0;
                if (flow > 0.0)
                {
                    water_a = (QTheta[layer, num_Qpoints] - hydraulicModel.get_theta(layer, psiad)) * QWeight[layer, num_Qpoints] * soilPhysical.Thickness[layer];
                    // TODO: no need to limit flow to water_d.
                    // water_d = (soilPhysical.SAT[layer + 1] - QTheta[layer + 1, 0]) * QWeight[layer + 1, 0] * soilPhysical.Thickness[layer + 1];
                    max_flow = soilPhysical.KS[layer] + soilPhysical.KS[layer + 1];
                    // max_flow = Math.Min(max_flow, Math.Min(water_a, water_d));
                    max_flow = Math.Min(max_flow, water_a);
                    max_flow = Math.Min(max_flow, flow);
                    flow = max_flow;
                }
                else
                {
                    //water_a = (QTheta[layer + 1, 0] - SimpleTheta(layer + 1, psiad)) * QWeight[layer + 1, 0] * soilPhysical.Thickness[layer + 1];
                    water_a = (QTheta[layer + 1, 1] - hydraulicModel.get_theta(layer + 1, psiad)) * QWeight[layer + 1, 1] * soilPhysical.Thickness[layer + 1];
                    water_d = (soilPhysical.SAT[layer] - QTheta[layer, num_Qpoints]) * QWeight[layer, num_Qpoints] * soilPhysical.Thickness[layer];
                    max_flow = soilPhysical.KS[layer];
                    max_flow = -Math.Min(max_flow, Math.Min(water_a, water_d));
                    max_flow = Math.Max(max_flow, flow);
                    flow = max_flow;
                }

                old_flow = flow;
                for (int iteration = 0; iteration < 100; ++iteration)
                {
                    theta_1 = QTheta[layer, num_Qpoints] - flow / (QWeight[layer, num_Qpoints] * soilPhysical.Thickness[layer]);
                    theta_2 = QTheta[layer + 1, 1] + flow / (QWeight[layer + 1, 1] * soilPhysical.Thickness[layer + 1]);
                    psi_1 = hydraulicModel.get_h(layer, theta_1);
                    psi_2 = hydraulicModel.get_h(layer + 1, theta_2);
                    grad = (psi_1 - psi_2) / depth + 1.0;

                    if (grad * flow > 0)
                    {
                        if (iteration == 0)
                            return flow;
                        else
                        {
                            min_flow = flow;
                            flow = (flow + max_flow) / 2;
                        }
                    }
                    else
                    {
                        max_flow = flow;
                        flow = (flow + min_flow) / 2;
                    }
                    double delta = Math.Abs(flow - old_flow);
                    if (delta < tolerance)
                        return flow;
                    old_flow = flow;
                }

                return flow;
            }
            else
            {
                double flow_d = 1.0;
                double flow_cum = 0.0;
                double theta;
                double psi;
                double K;

                theta = QTheta[layer, num_Qpoints];
                flow = QK[layer, num_Qpoints];
                //water_a = (QTheta[layer, 2] - SimpleTheta(layer, psiad)) * QWeight[layer, 2] * soilPhysical.Thickness[layer];
                water_a = (QTheta[layer, num_Qpoints] - hydraulicModel.get_theta(layer, psiad)) * QWeight[layer, num_Qpoints] * soilPhysical.Thickness[layer];
                flow = Math.Min(flow, water_a);

                if (flow <= flow_d)
                {
                    return flow;
                }

                for (int i = 0; i < 100; ++i)
                {
                    if (flow <= flow_d)
                    {
                        flow_cum += flow;
                        if (flow_cum >= QK[layer, num_Qpoints])
                            return QK[layer, num_Qpoints];
                        else
                            return flow_cum;
                    }

                    flow = flow_d;
                    flow_cum += flow;
                    if (flow_cum >= QK[layer, num_Qpoints])
                    {
                        return QK[layer, num_Qpoints];
                    }
                    
                    theta -= flow / (QWeight[layer, num_Qpoints] * soilPhysical.Thickness[layer]);
                    //psi = Suction(layer, theta);
                    //K = SimpleK(layer, psi);
                    psi = hydraulicModel.get_h(layer, theta);
                    K = hydraulicModel.get_K(layer, psi);

                    flow = K;
                }

                return flow_cum;
            }
        }

        /// <summary>
        /// Reset status for quadrature points to old values.
        /// </summary>
        /// <param name="layer"></param>
        private void ResetQ(int layer)
        {

        }

        /// <summary>
        /// Calculate the available water for gravitational flow from this layer.
        /// </summary>
        /// <param name="layer"></param>
        /// <returns></returns>
        private double AvailableWater(int layer)
        {
            double availableWater = 0.0;

            for (int q = 0; q < num_Qpoints; ++q)
            {
                availableWater += (QTheta[layer, q] - soilPhysical.DUL[layer]) * QWeight[layer, q] * soilPhysical.Thickness[layer];
            }

            return Math.Max(0.0, availableWater);
        }
    }
}
