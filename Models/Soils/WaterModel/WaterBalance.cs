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

        // Constants
        const double effpar = 0.184;

        // Units in mm
        const double psi_ll15 = -150000.0;
        const double psiad = -1e7;
        const double psi0 = -0.6e8;

        /// <summary>Number of soil layers.</summary>
        private int num_layers { get { return soilPhysical.Thickness.Length; } }
        private double[] psid;

        // Intermediate variables for soil hydraulic curves
        private double[,] DELk;
        private double[,] Mk;
        private double[,] M0;
        private double[,] M1;
        private double[,] Y0;
        private double[,] Y1;
        private double[] MicroP;
        private double[] MicroKs;
        private double[] Kdula;
        private double[] kdul;
        private double[] MacroP;

        // Parameters at quadrature points
        // Consider encapsulate all parameters in a class/struct
        private string QuadratureRule;
        private int num_Qpoints;
        private double[,] QDepth;
        private double[,] QWeight;
        private double[,] QTheta;
        private double[,] Qpsi;
        private double[,] QK;

        private double[] OldWater;
        private double[] InterfaceFlow;
        private int[] FlowType;


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

        // --- Event handlers ------------------------------------------------------------

        /// <summary>Called when a simulation commences.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The event data.</param>
        [EventSubscribe("Commencing")]
        private void OnStartOfSimulation(object sender, EventArgs e)
        {
            Initialise();
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
                    //Water[irrigationLayer] += irrigation.Amount;
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
            // double es = evaporationModel.Calculate();
            // Water[0] = Water[0] - es;

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


        /// Interpolate hydraulic parameters from APSIMX inputs
        /// Code based on SWIM
        /// TODO: fix layer and node

        ///<summary>Perform initial calculations for hydraulic curves</summary>
        private void InitCalc()
        {
            //+  Purpose
            //   Perform initial calculations from input parameters and prepare for simulation

            // change units of params to normal SWIM units
            // ie. cm and hours etc.
            //InitChangeUnits();

            // ------------------- CALCULATE CURRENT TIME -------------------------
            //int time_mins = TimeToMins(apsim_time);
            //t = Time(year, day, time_mins);

            // ----------------- SET UP NODE SPECIFICATIONS -----------------------

            // safer to use number returned from read routine
            //int num_layers = soilPhysical.Thickness.Length;
            //if (n != num_layers - 1)
            //    ResizeProfileArrays(num_layers);

            //for (int i = 0; i <= n; i++)
            //    dx[i] = soilPhysical.Thickness[i] / 10.0;

            //x[0] = 0.0;
            //x[1] = 2.0 * dx[0] + x[0];

            //for (int i = 1; i < n; i++)
            //    x[i] = MathUtilities.Sum(dx, 0, i - 1) + dx[i] / 2;

            //x[n] = MathUtilities.Sum(dx);

            //      p%dx(0) = 0.5*(p%x(1) - p%x(0))
            //      do 10 i=1,p%n-1
            //         p%dx(i) = 0.5*(p%x(i+1)-p%x(i-1))
            //   10 continue
            //      p%dx(p%n) = 0.5*(p%x(p%n)-p%x(p%n-1))


            // ------- IF USING SIMPLE SOIL SPECIFICATION CALCULATE PROPERTIES -----

            DELk = new double[num_layers, 4];
            Mk = new double[num_layers, 4];
            M0 = new double[num_layers, 5];
            M1 = new double[num_layers, 5];
            Y0 = new double[num_layers, 5];
            Y1 = new double[num_layers, 5];
            MicroP = new double[num_layers];
            MicroKs = new double[num_layers];
            kdul = new double[num_layers];
            Kdula = new double[num_layers];
            MacroP = new double[num_layers];
            psid = new double[num_layers];

            InterfaceFlow = new double[num_layers + 1];
            FlowType = new int[num_layers + 1];
            OldWater = new double[num_layers];

            OldWater = Water;

            // TODO: set KDul as constant and PSIDul as variable for different soils
            PSIDul = -3400.0;
            KDul = 0.1;

            for (int layer = 0; layer < num_layers; ++layer)
            {
                kdul[layer] = KDul;
                psid[layer] = PSIDul;
            }

            SetupThetaCurve();
            SetupKCurve();

            QuadratureRule = "Gaussian";
            num_Qpoints = 3;

            InitQuadrature(QuadratureRule, num_Qpoints);       
        }

        private void SetupThetaCurve()
        {
            for (int layer = 0; layer < num_layers; layer++)
            {
                DELk[layer, 0] = (soilPhysical.DUL[layer] - soilPhysical.SAT[layer]) / (Math.Log10(-psid[layer]));
                DELk[layer, 1] = (soilPhysical.LL15[layer] - soilPhysical.DUL[layer]) / (Math.Log10(-psi_ll15) - Math.Log10(-psid[layer]));
                DELk[layer, 2] = -soilPhysical.LL15[layer] / (Math.Log10(-psi0) - Math.Log10(-psi_ll15));
                DELk[layer, 3] = -soilPhysical.LL15[layer] / (Math.Log10(-psi0) - Math.Log10(-psi_ll15));

                Mk[layer, 0] = 0.0;
                Mk[layer, 1] = (DELk[layer, 0] + DELk[layer, 1]) / 2.0;
                Mk[layer, 2] = (DELk[layer, 1] + DELk[layer, 2]) / 2.0;
                Mk[layer, 3] = DELk[layer, 3];

                // First bit might not be monotonic so check and adjust
                double alpha = Mk[layer, 0] / DELk[layer, 0];
                double beta = Mk[layer, 1] / DELk[layer, 0];
                double phi = alpha - (Math.Pow(2.0 * alpha + beta - 3.0, 2.0) / (3.0 * (alpha + beta - 2.0)));
                if (phi <= 0)
                {
                    double tau = 3.0 / Math.Sqrt(alpha * alpha + beta * beta);
                    Mk[layer, 0] = tau * alpha * DELk[layer, 0];
                    Mk[layer, 1] = tau * beta * DELk[layer, 0];
                }

                M0[layer, 0] = 0.0;
                M1[layer, 0] = 0.0;
                Y0[layer, 0] = soilPhysical.SAT[layer];
                Y1[layer, 0] = soilPhysical.SAT[layer];

                M0[layer, 1] = Mk[layer, 0] * (Math.Log10(-psid[layer]) - 0.0);
                M1[layer, 1] = Mk[layer, 1] * (Math.Log10(-psid[layer]) - 0.0);
                Y0[layer, 1] = soilPhysical.SAT[layer];
                Y1[layer, 1] = soilPhysical.DUL[layer];

                M0[layer, 2] = Mk[layer, 1] * (Math.Log10(-psi_ll15) - Math.Log10(-psid[layer]));
                M1[layer, 2] = Mk[layer, 2] * (Math.Log10(-psi_ll15) - Math.Log10(-psid[layer]));
                Y0[layer, 2] = soilPhysical.DUL[layer];
                Y1[layer, 2] = soilPhysical.LL15[layer];

                M0[layer, 3] = Mk[layer, 2] * (Math.Log10(-psi0) - Math.Log10(-psi_ll15));
                M1[layer, 3] = Mk[layer, 3] * (Math.Log10(-psi0) - Math.Log10(-psi_ll15));
                Y0[layer, 3] = soilPhysical.LL15[layer];
                Y1[layer, 3] = 0.0;

                M0[layer, 4] = 0.0;
                M1[layer, 4] = 0.0;
                Y0[layer, 4] = 0.0;
                Y1[layer, 4] = 0.0;
            }
        }

        private void SetupKCurve()
        {
            for (int layer = 0; layer < num_layers; layer++)
            {
                double b = -Math.Log(PSIDul / psi_ll15) / Math.Log(soilPhysical.DUL[layer] / soilPhysical.LL15[layer]);
                MicroP[layer] = b * 2.0 + 3.0;
                Kdula[layer] = Math.Min(0.99 * kdul[layer], soilPhysical.KS[layer]);
                MicroKs[layer] = Kdula[layer] / Math.Pow(soilPhysical.DUL[layer] / soilPhysical.SAT[layer], MicroP[layer]);

                double Sdul = soilPhysical.DUL[layer] / soilPhysical.SAT[layer];
                MacroP[layer] = Math.Log10(Kdula[layer] / 99.0 / (soilPhysical.KS[layer] - MicroKs[layer])) / Math.Log10(Sdul);
            }
        }

        private double Suction(int node, double theta)
        {
            //  Purpose
            //  Calculate the suction for a given water content for a given node.
            // Calculation from SWIM3 is wrong.
            // TODO: temporary fix; need further thinking.
            const int maxIterations = 1000;
            const double tolerance = 1e-9;
            // const double dpsi = 0.01;
            double dpsi = 0.01;

            double theta_up;
            double theta_low;
            double psi_up;
            double psi_low;
            double delta_old = dpsi;

            // Temporary fix.
            psi_up = 0.0;
            theta_up = soilPhysical.SAT[node];
            psi_low = psid[node];
            theta_low = SimpleTheta(node, psi_low);
            if (theta < theta_low)
            {
                psi_up = psi_low;
                theta_up = theta_low;
                psi_low = psi_ll15;
                theta_low = SimpleTheta(node, psi_low);
                if (theta < theta_low)
                {
                    psi_up = psi_low;
                    theta_up = theta_low;
                    psi_low = psi0;
                    theta_low = SimpleTheta(node, psi_low);
                }
            }

            if (theta >= soilPhysical.SAT[node])
                return 0.0;
            else
            {
                double psiValue = -3400.0; 
                for (int iter = 0; iter < maxIterations; iter++)
                {
                    double est = SimpleTheta(node, psiValue);
                    double m = (SimpleTheta(node, psiValue + dpsi) - est) / dpsi;
                    double delta = (est - theta) / Math.Abs(est - theta) * Math.Max(dpsi, Math.Abs((est - theta) / m));
                    if (MathUtilities.FloatsAreEqual(delta, -delta_old))
                    {
                        delta /= 2.0;
                        dpsi /= 2.0;
                    }
                    delta_old = delta;
                    double psi_temp = psiValue - delta;

                    if (Math.Abs(est - theta) < tolerance)
                        break;
                    // psiValue -= Math.Min(-dpsi, (est - theta) / m);
                    if (psi_temp > psi_up | psi_temp < psi_low)
                    {
                        if (est < theta)
                        {
                            psi_low = psiValue;
                            theta_low = est;
                            psiValue = (psiValue + psi_up) / 2;
                        }
                        else
                        {
                            psi_up = psiValue;
                            theta_up = est;
                            psiValue = (psiValue + psi_low) / 2;
                        }
                    }
                    else
                        psiValue = psi_temp;
                }
                if (psiValue > 0.0)
                    System.Console.WriteLine("Error in psi.");
                return psiValue;
            }
        }

        private double SimpleS(int layer, double psiValue)
        {
            //  Purpose
            //      Calculate S for a given node for a specified suction.
            return SimpleTheta(layer, psiValue) / soilPhysical.SAT[layer];
        }

        private double SimpleTheta(int layer, double psiValue)
        {
            //  Purpose
            //     Calculate Theta for a given node for a specified suction.
            int i;
            double t;

            if (psiValue >= -1.0)
            {
                i = 0;
                t = 0.0;
            }
            else if (psiValue > psid[layer])
            {
                i = 1;
                t = (Math.Log10(-psiValue) - 0.0) / (Math.Log10(-psid[layer]) - 0.0);
            }
            else if (psiValue > psi_ll15)
            {
                i = 2;
                t = (Math.Log10(-psiValue) - Math.Log10(-psid[layer])) / (Math.Log10(-psi_ll15) - Math.Log10(-psid[layer]));
            }
            else if (psiValue > psi0)
            {
                i = 3;
                t = (Math.Log10(-psiValue) - Math.Log10(-psi_ll15)) / (Math.Log10(-psi0) - Math.Log10(-psi_ll15));
            }
            else
            {
                i = 4;
                t = 0.0;
            }

            double tSqr = t * t;
            double tCube = tSqr * t;

            return (2 * tCube - 3 * tSqr + 1) * Y0[layer, i] + (tCube - 2 * tSqr + t) * M0[layer, i]
                    + (-2 * tCube + 3 * tSqr) * Y1[layer, i] + (tCube - tSqr) * M1[layer, i];
        }

        private void Interp(int node, double tpsi, out double tth, out double thd, out double hklg, out double hklgd)
        {
            //  Purpose
            //   interpolate water characteristics for given potential for a given
            //   node.

            const double dpsi = 0.0001;
            double temp;

            tth = SimpleTheta(node, tpsi);
            temp = SimpleTheta(node, tpsi + dpsi);
            thd = (temp - tth) / Math.Log10((tpsi + dpsi) / tpsi);
            hklg = Math.Log10(SimpleK(node, tpsi));
            temp = Math.Log10(SimpleK(node, tpsi + dpsi));
            hklgd = (temp - hklg) / Math.Log10((tpsi + dpsi) / tpsi);
        }

        private double SimpleK(int layer, double psiValue)
        {
            //  Purpose
            //      Calculate Conductivity for a given node for a specified suction.

            double S = SimpleS(layer, psiValue);
            double simpleK;

            if (S <= 0.0)
                simpleK = 1e-100;
            else
            {
                double microK = MicroKs[layer] * Math.Pow(S, MicroP[layer]);

                if (MicroKs[layer] >= soilPhysical.KS[layer])
                    simpleK = microK;
                else
                {
                    double macroK = (soilPhysical.KS[layer] - MicroKs[layer]) * Math.Pow(S, MacroP[layer]);
                    simpleK = microK + macroK;
                }
            }
            return simpleK;
            //return simpleK / 24.0 / 10.0;
        }

        private double Theta(int node, double suction)
        {
            double theta;
            double thd;
            double hklg;
            double hklgd;

            Interp(node, suction, out theta, out thd, out hklg, out hklgd);
            return theta;
        }

        /// <summary>
        /// Setup quadrature points for each layer. Default is 3 points Gaussian quadrature.
        /// </summary>
        private void InitQuadrature(string quadrature_rule = "Gaussian", int num_points = 3)
        {
            double[] y;
            double[] weight;
            y = new double[num_points];
            weight = new double[num_points];

            QDepth = new double[num_layers, num_points];
            QWeight = new double[num_layers, num_points];
            QTheta = new double[num_layers, num_points];
            QK = new double[num_layers, num_points];
            Qpsi = new double[num_layers, num_points];

            if (quadrature_rule == "Gaussian")
            {
                if (num_points == 3)
                {
                    y[0] = (1 - Math.Sqrt(3.0 / 5.0)) / 2.0;
                    y[1] = 1.0 / 2.0;
                    y[2] = (1 + Math.Sqrt(3.0 / 5.0)) / 2.0;
                    weight[0] = 5.0 / 9.0 / 2.0;
                    weight[1] = 8.0 / 9.0 / 2.0;
                    weight[2] = 5.0 / 9.0 / 2.0;

                    for (int layer = 0; layer < num_layers; ++layer)
                    {
                        for (int point = 0; point < num_points; ++point)
                        {
                            QTheta[layer, point] = SW[layer];
                            Qpsi[layer, point] = Suction(layer, QTheta[layer, point]);
                            QK[layer, point] = SimpleK(layer, Qpsi[layer, point]);
                            QDepth[layer, point] = soilPhysical.ThicknessCumulative[layer] + (y[point] - 1) * soilPhysical.Thickness[layer];
                            QWeight[layer, point] = weight[point];
                        }
                    }
                }
            }
        }

        /// <summary> Calculate water flow at each interface.</summary>
        private void WaterFlow()
        {
            // TODO: include Pond calculation.
            // Check flow types at each interfaces
            double[] NewWater = new double[num_layers];
            double[] Source = new double[num_layers];
            double[] BackFlow = new double[num_layers];
            double[] ExcessW = new double[num_layers];
            bool[] is_Saturated = new bool[num_layers];
            double MaxFlow;
            double Theta;
            double UpFlow;
            double DownFlow;

            for (int layer = 0; layer < num_layers; ++layer)
            {
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
            FlowType[0] = 1;

            // Check for sources in the profile, e.g. subsurface irrigation
            // TODO: need to consider subsurface irrigation in the first layer. Use -1 for surface irrigation and 
            foreach (var irrigation in irrigations)
            {
                if (irrigation.Amount > 0)
                {
                    int irrigationLayer = SoilUtilities.LayerIndexOfDepth(soilPhysical.Thickness, Convert.ToInt32(irrigation.Depth, CultureInfo.InvariantCulture));
                    if (irrigationLayer == 0)
                        Infiltration += irrigation.Amount;
                    else
                    {
                        Source[irrigationLayer] += irrigation.Amount;
                    }
                }
            }

            // TODO: compare infiltration with potential evaporation and subtract one from another

            // TODO: drainage boundary assumed for the bottom; add other possibilities.

            for (int layer = 0; layer < num_layers; ++layer)
            {
                NewWater[layer] = Water[layer];
                if (Source[layer] > 0.0)
                {
                    NewWater[layer] += Source[layer];
                    // TODO: may just use ExcessW to indicate saturation
                    is_Saturated[layer] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer], soilPhysical.SATmm[layer]);
                    ExcessW[layer] = Math.Max(0.0, NewWater[layer] - soilPhysical.SATmm[layer]);
                    Redistribute(layer, Source[layer], 0);
                }
            }

            if (Infiltration > 0.0)
            {
                // TODO: This flux might limit infiltration when soil is relatively dry.
                DownFlow = QK[0, 0] * ((Pond - (Qpsi[0, 0]) / QDepth[0, 0] + 1));
                MaxFlow = Math.Max(DownFlow, 2 * soilPhysical.KS[0]);
                InterfaceFlow[0] = Math.Min(MaxFlow, Infiltration);
                NewWater[0] += InterfaceFlow[0];
                Redistribute(0, InterfaceFlow[0], -1);
                if (ExcessW[0] > 0.0)
                    ExcessW[0] += InterfaceFlow[0];
                else
                {
                    is_Saturated[0] = MathUtilities.IsGreaterThanOrEqual(NewWater[0], soilPhysical.SATmm[0]);
                    ExcessW[0] = Math.Max(0.0, NewWater[0] - soilPhysical.SATmm[0]);
                }

            }

            //for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    NewWater[layer] = Water[layer] + Source[layer];
            //    // TODO: may just use ExcessW to indicate saturation
            //    is_Saturated[layer] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer], soilPhysical.SATmm[layer]);
            //    ExcessW[layer] = Math.Max(0.0, NewWater[layer] - soilPhysical.SATmm[layer]);
            //    Redistribute(layer);
            //}

            for (int layer = 0; layer < num_layers - 1; ++layer)
            {
                if (FlowType[layer] == 1)
                {
                    //if (ExcessW[layer] > 0.0)
                    //    ExcessW[layer] += InterfaceFlow[layer];
                    //else
                    //{
                    //    NewWater[layer] += InterfaceFlow[layer];
                    //    is_Saturated[layer] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer], soilPhysical.SATmm[layer]);
                    //    ExcessW[layer] = Math.Max(0.0, NewWater[layer] - soilPhysical.SATmm[layer]);
                    //}

                    if (ExcessW[layer] > 0.0)
                    {
                        FlowType[layer + 1] = 1;
                        InterfaceFlow[layer + 1] = ExcessW[layer];
                        MaxFlow = 2 * soilPhysical.KS[layer + 1];
                        BackFlow[layer] = InterfaceFlow[layer + 1] - MaxFlow;
                        if (BackFlow[layer] > 0.0)
                        {
                            InterfaceFlow[layer + 1] -= BackFlow[layer];
                            Redistribute(layer + 1, InterfaceFlow[layer + 1], -1);
                        }
                        else
                        {
                            Redistribute(layer + 1, InterfaceFlow[layer + 1], -1);
                            DownFlow = UnsatFlow(layer);
                            if (MathUtilities.IsGreaterThanOrEqual(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
                                DownFlow = MaxFlow - InterfaceFlow[layer + 1];

                            InterfaceFlow[layer + 1] += DownFlow;
                            Redistribute(layer, DownFlow, 1);
                            Redistribute(layer + 1, DownFlow, -1);
                        }
                        NewWater[layer] -= InterfaceFlow[layer + 1];
                        NewWater[layer + 1] += InterfaceFlow[layer + 1];
                        if (ExcessW[layer + 1] > 0.0)
                            ExcessW[layer + 1] += InterfaceFlow[layer + 1];
                        else
                        {
                            is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
                            ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                        }
                    }
                    else
                    {
                        DownFlow = UnsatFlow(layer);
                        InterfaceFlow[layer + 1] = DownFlow;
                        if (DownFlow > 0.0)
                        {
                            FlowType[layer + 1] = 1;
                            Redistribute(layer, DownFlow, 1);
                            NewWater[layer] -= InterfaceFlow[layer + 1];
                            Redistribute(layer + 1, DownFlow, -1);
                            NewWater[layer + 1] += InterfaceFlow[layer + 1];
                            if (ExcessW[layer + 1] > 0.0)
                                ExcessW[layer + 1] += InterfaceFlow[layer + 1];
                            else
                            {
                                is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
                                ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                            }
                        }
                    }
                }
                else
                {
                    if (ExcessW[layer] > 0.0)
                    {
                        double amount;
                        bool is_upflow = true;
                        FlowType[layer + 1] = 1;
                        InterfaceFlow[layer] = 0.0;
                        InterfaceFlow[layer + 1] = 0.0;

                        do
                        {
                            if (ExcessW[layer] < 2.0)
                                amount = ExcessW[layer];
                            else
                                amount = 2.0;
                            ExcessW[layer] -= 2.0;

                            if (is_upflow)
                            {
                                // This layer is fully saturated. The downward flow is calculate assuming pressure building in the whole layer.
                                double up = -QK[layer - 1, 2] * ((Qpsi[layer - 1, 2] - 0.0) / (soilPhysical.ThicknessCumulative[layer - 1] - QDepth[layer - 1, 2]) + 1);
                                double down = QK[layer + 1, 0] * ((0.0 - Qpsi[layer + 1, 0]) / (QDepth[layer + 1, 0] - soilPhysical.ThicknessCumulative[layer]) + 2);
                                if (up > 0)
                                {
                                    UpFlow = up / (up + down) * amount;
                                    // Redistribute(layer - 1, -UpFlow, -1);
                                    Theta = QTheta[layer - 1, 2] + UpFlow / (soilPhysical.Thickness[layer - 1] * QWeight[layer - 1, 2]);
                                    if (MathUtilities.IsGreaterThanOrEqual(Theta, soilPhysical.SAT[layer - 1]))
                                    {
                                        UpFlow -= (Theta - soilPhysical.SAT[layer - 1]) * (soilPhysical.Thickness[layer - 1] * QWeight[layer - 1, 2]);
                                        is_upflow = false;
                                        // Theta = soilPhysical.SAT[layer - 1];
                                    }
                                    InterfaceFlow[layer] += -UpFlow;
                                    Redistribute(layer - 1, -UpFlow, 1);
                                    // QTheta[layer - 1, 2] = Theta;
                                    Qpsi[layer - 1, 2] = Suction(layer - 1, QTheta[layer - 1, 2]);
                                    QK[layer - 1, 2] = SimpleK(layer - 1, Qpsi[layer - 1, 2]);
                                }
                                else
                                {
                                    is_upflow = false;
                                    UpFlow = 0.0;
                                }
                            }
                            else
                                UpFlow = 0.0;

                            DownFlow = amount - UpFlow;
                            //if (!is_sat_down)
                            //{
                            //    Theta = QTheta[layer + 1, 0] + DownFlow / (soilPhysical.Thickness[layer + 1] * QWeight[layer + 1, 0]);
                            //    if (MathUtilities.IsGreaterThan(Theta, soilPhysical.SAT[layer + 1]))
                            //    {
                            //        Theta = soilPhysical.SAT[layer + 1];
                            //        is_sat_down = true;
                            //    }
                            //    Redistribute(layer + 1, DownFlow, 1);
                            //    // QTheta[layer + 1, 0] = Theta;
                            //    Qpsi[layer + 1, 0] = Suction(layer + 1, QTheta[layer + 1, 0]);
                            //    QK[layer + 1, 0] = SimpleK(layer + 1, Qpsi[layer + 1, 0]);
                            //}

                            InterfaceFlow[layer + 1] += DownFlow;
                            Redistribute(layer + 1, DownFlow, -1);
                        } while (ExcessW[layer] > 0.0);

                        // continue unsatflow
                        if (is_upflow)
                        {
                            UpFlow = -UnsatFlow(layer - 1);
                            if (UpFlow > 0.0)
                            {
                                InterfaceFlow[layer] += -UpFlow;
                                Redistribute(layer - 1, -UpFlow, 1);
                                Redistribute(layer, -UpFlow, -1);
                            }
                        }

                        MaxFlow = 2 * soilPhysical.KS[layer + 1];
                        if (MathUtilities.IsLessThan(InterfaceFlow[layer + 1], MaxFlow))
                        {
                            DownFlow = UnsatFlow(layer);
                            if (MathUtilities.IsGreaterThan(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
                                DownFlow = MaxFlow - InterfaceFlow[layer + 1];
                            InterfaceFlow[layer + 1] += DownFlow;
                            Redistribute(layer, DownFlow, 1);
                        }

                        NewWater[layer - 1] -= InterfaceFlow[layer];
                        NewWater[layer] += InterfaceFlow[layer];
                        NewWater[layer] -= InterfaceFlow[layer + 1];
                        NewWater[layer + 1] += InterfaceFlow[layer + 1];
                        if (ExcessW[layer + 1] > 0.0)
                            ExcessW[layer + 1] += InterfaceFlow[layer + 1];
                        else
                        {
                            is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
                            ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                        }
                    }
                    else
                    {
                        UpFlow = -UnsatFlow(layer - 1);
                        InterfaceFlow[layer] = -UpFlow;
                        Redistribute(layer - 1, -UpFlow, 1);
                        Redistribute(layer, -UpFlow, -1);
                        NewWater[layer - 1] -= InterfaceFlow[layer];
                        NewWater[layer] += InterfaceFlow[layer];

                        DownFlow = UnsatFlow(layer);
                        InterfaceFlow[layer + 1] = DownFlow;
                        if (DownFlow > 0.0)
                        {
                            FlowType[layer + 1] = 1;
                            Redistribute(layer, DownFlow, 1);
                            Redistribute(layer + 1, DownFlow, -1);
                            NewWater[layer] -= InterfaceFlow[layer + 1];
                            NewWater[layer + 1] += InterfaceFlow[layer + 1];
                            if (ExcessW[layer + 1] > 0.0)
                                ExcessW[layer + 1] += InterfaceFlow[layer + 1];
                            else
                            {
                                is_Saturated[layer + 1] = MathUtilities.IsGreaterThanOrEqual(NewWater[layer + 1], soilPhysical.SATmm[layer + 1]);
                                ExcessW[layer + 1] = Math.Max(0.0, NewWater[layer + 1] - soilPhysical.SATmm[layer + 1]);
                            }
                        }
                    }
                }
            }

            // Last layer.
            // TODO: other boundary type.
            {
                int layer = num_layers- 1;
                InterfaceFlow[layer + 1] = 0.0;

                if (FlowType[layer] == 1)
                {
                    if (ExcessW[layer] > 0.0)
                    {
                        InterfaceFlow[layer + 1] = ExcessW[layer];
                        MaxFlow = soilPhysical.KS[layer];
                        BackFlow[layer] = InterfaceFlow[layer + 1] - MaxFlow;
                        if (BackFlow[layer] > 0.0)
                        {
                            InterfaceFlow[layer + 1] -= BackFlow[layer];
                        }
                        else
                        {
                            DownFlow = UnsatFlow(layer);
                            if (MathUtilities.IsGreaterThan(InterfaceFlow[layer + 1] + DownFlow, MaxFlow))
                                DownFlow = MaxFlow - InterfaceFlow[layer + 1];
                            InterfaceFlow[layer + 1] += DownFlow;
                            Redistribute(layer, DownFlow, 1);
                        }
                    }
                    else
                    {
                        DownFlow = UnsatFlow(layer);
                        InterfaceFlow[layer + 1] = DownFlow;
                        Redistribute(layer, QK[layer, 2], 1);
                    }
                }
                else
                {
                    UpFlow = -UnsatFlow(layer - 1);
                    InterfaceFlow[layer] = -UpFlow;
                    Redistribute(layer - 1, -UpFlow, 1);
                    Redistribute(layer, -UpFlow, -1);
                    NewWater[layer - 1] -= InterfaceFlow[layer];
                    NewWater[layer] += InterfaceFlow[layer];

                    DownFlow = UnsatFlow(layer);
                    InterfaceFlow[layer + 1] = DownFlow;
                    Redistribute(layer, DownFlow, 1);
                }

                NewWater[layer] -= InterfaceFlow[layer + 1];
            }

            // Backup flow
            for (int layer = num_layers - 1; layer > 0; --layer)
            {
                if (MathUtilities.IsGreaterThan(NewWater[layer], soilPhysical.SATmm[layer]))
                {
                    BackFlow[layer] = NewWater[layer] - soilPhysical.SATmm[layer];
                    InterfaceFlow[layer] -= BackFlow[layer];
                    Redistribute(layer - 1, -BackFlow[layer], 1);
                    NewWater[layer - 1] += BackFlow[layer];
                }
            }
            if (MathUtilities.IsGreaterThan(NewWater[0], soilPhysical.SATmm[0]))
            {
                BackFlow[0] = NewWater[0] - soilPhysical.SATmm[0];
                InterfaceFlow[0] -= BackFlow[0];
            }

            // Unsaturated flow within each layer.
            for (int layer = 0; layer < num_layers; ++layer)
            {
                Redistribute(layer, 1);
            }

            //// Previous method ---------------------------------------------------------------------------------------------------------
            //    for (int layer = 0; layer < num_layers; ++layer)
            //{
            //    if (layer == num_layers -1)
            //    {
            //        // Free drainage boundary
            //        InterfaceFlow[layer + 1] = (QK[layer, 2]);
            //    }
            //    else
            //    {
            //        InterfaceFlow[layer + 1] = (QK[layer, 2] + QK[layer + 1, 0]) / 2
            //                                 * ((Qpsi[layer, 2] - Qpsi[layer + 1, 0]) / (QDepth[layer + 1, 0] - QDepth[layer, 2]) + 1);
            //        // MaxFlow = (QK[layer, 2] + QK[layer + 1, 0]);
            //        MaxFlow = soilPhysical.KS[layer] + soilPhysical.KS[layer + 1];
            //        if (InterfaceFlow[layer + 1] > 0)
            //        {
            //            double MaxDownSat = Water[layer] + Source[layer] + InterfaceFlow[layer] - soilPhysical.DULmm[layer];
            //            MaxDownSat = Math.Max(0.0, MaxDownSat);
            //            // TODO: update psi at lower layer quadrature point
            //            double MaxDownUnsat = kdul[layer] * ((psid[layer] - Qpsi[layer + 1, 0]) / (QDepth[layer + 1, 0] - QDepth[layer, 2]) + 1);
            //            MaxDownUnsat = Math.Max(0.0, MaxDownUnsat);
            //            double MaxDownFlow = MaxDownSat + MaxDownUnsat;
            //            MaxFlow = Math.Min(MaxFlow, MaxDownFlow);
            //            InterfaceFlow[layer + 1] = Math.Min(InterfaceFlow[layer + 1], MaxFlow);
            //        }
            //        else
            //        {
            //            double MaxUpSat = -(Water[layer + 1] + Source[layer + 1] - soilPhysical.DULmm[layer + 1]);
            //            MaxUpSat = Math.Min(0.0, MaxUpSat);
            //            double MaxUpUnsat = kdul[layer + 1] * ((Qpsi[layer, 2] - psid[layer + 1]) / (QDepth[layer + 1, 0] - QDepth[layer, 2]) + 1);
            //            MaxUpUnsat = Math.Min(0.0, MaxUpUnsat);
            //            double MaxUpFlow = MaxUpSat + MaxUpUnsat;
            //            MaxFlow = Math.Max(-MaxFlow, MaxUpFlow);
            //            InterfaceFlow[layer + 1] = Math.Max(InterfaceFlow[layer + 1], MaxFlow);
            //        }
            //    }
                
            //    // TODO: make sure enough water is available to flow out from the region.
            //    ExcessW[layer] = Water[layer] + Source[layer] + InterfaceFlow[layer] - InterfaceFlow[layer + 1] - soilPhysical.SATmm[layer];

            //    if(ExcessW[layer] > 0.0)
            //    {
            //        is_Saturated[layer] = true;
            //        if (layer != 0 & FlowType[layer] == 0)
            //        {
            //            InterfaceFlow[layer] = (QK[layer - 1, 2])
            //                                   * ((Qpsi[layer - 1, 2] - 0.0) / (soilPhysical.ThicknessCumulative[layer - 1] - QDepth[layer - 1, 2]) + 1);
            //            // TODO: *** need to limit this flow; or delete this step; if source[layer] is 0, no need to do this step; or set FlowType to 1 for scenarios
            //        }
            //        ExcessW[layer] = Water[layer] + Source[layer] + InterfaceFlow[layer] - InterfaceFlow[layer + 1] - soilPhysical.SATmm[layer];
            //        if (layer != num_layers - 1)
            //            MaxFlow = 2 * soilPhysical.KS[layer + 1];
            //        else
            //            MaxFlow = 2 * soilPhysical.KS[layer];

            //        if (ExcessW[layer] > MaxFlow)
            //        {
            //            InterfaceFlow[layer + 1] = MaxFlow;
            //            BackFlow[layer] = ExcessW[layer] - MaxFlow;
            //        }
            //        else if (ExcessW[layer] > InterfaceFlow[layer + 1])
            //        {
            //            InterfaceFlow[layer + 1] = ExcessW[layer];
            //        }
            //    }
            //}

            //for (int layer = num_layers - 1; layer >= 0; --layer)
            //{
            //    if (is_Saturated[layer])
            //    {
            //        if (layer == num_layers - 1)
            //        {
            //            BackFlow[layer] = 0;
            //        }
            //        else
            //        {
            //            BackFlow[layer] = Water[layer] + Source[layer] + InterfaceFlow[layer] - InterfaceFlow[layer + 1] - soilPhysical.SATmm[layer];
            //        }
                    
            //        InterfaceFlow[layer] -= BackFlow[layer];
            //    }
            //    else
            //    {
            //        // check if saturated
            //        ExcessW[layer] = Water[layer] + Source[layer] + InterfaceFlow[layer] - InterfaceFlow[layer + 1] - soilPhysical.SATmm[layer];
            //        if (ExcessW[layer] > 0)
            //        {
            //            is_Saturated[layer] = true;
            //            InterfaceFlow[layer] -= ExcessW[layer];
            //            if (InterfaceFlow[layer] < 0)
            //                BackFlow[layer] = -InterfaceFlow[layer];
            //        }
            //    }
            //}
        }

        /// <summary>
        /// Redistribute water within one layer.
        /// </summary>
        /// <param name="layer"></param>
        /// <param name="method"></param>
        private void Redistribute(int layer, int method = 0)
        {
            // This method is only used to calculate water distribution after changes outside of this class (root uptake).

            const double eff = 0.9;
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

                        for (int n = 0; n < num_Qpoints; ++n)
                        {
                            depth[n] = soilPhysical.Thickness[layer] * QWeight[layer, n];
                            QTheta[layer, n] -= uptake * QWeight[layer, n] / depth[n];
                        }

                        for (int n = 0; n < num_Qpoints; ++n)
                        {
                            Qpsi[layer, n] = Suction(layer, QTheta[layer, n]);
                            QK[layer, n] = SimpleK(layer, Qpsi[layer, n]);
                        }
                    }

                    break;
                case 1:
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
                    psi_mean = Suction(layer, theta_mean);
                    cap = (SimpleTheta(layer, psi_mean + delta) - SimpleTheta(layer, psi_mean - delta)) / (2 * delta);
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

                    for (int n = 0; n < num_Qpoints; ++n)
                    {
                        Qpsi[layer, n] = Suction(layer, QTheta[layer, n]);
                        QK[layer, n] = SimpleK(layer, Qpsi[layer, n]);
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
            switch (location)
            {
                case 1:
                    if (flow > 0.0)
                    {
                        // This flow shouldn't empty this section.
                        QTheta[layer, 2] -= flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsLessThan(QTheta[layer, 2], soilPhysical.AirDry[layer]))
                            System.Console.WriteLine("Too much water flow out from this section.");
                        Qpsi[layer, 2] = Suction(layer, QTheta[layer, 2]);
                        QK[layer, 2] = SimpleK(layer, Qpsi[layer, 2]);
                    }
                    else if (flow < 0.0)
                    {
                        QTheta[layer, 2] -= flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsGreaterThan(QTheta[layer, 2], soilPhysical.SAT[layer]))
                        {
                            flow = -(QTheta[layer, 2] - soilPhysical.SAT[layer]) * (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                            QTheta[layer, 2] = soilPhysical.SAT[layer];
                            Qpsi[layer, 2] = 0.0;
                            QK[layer, 2] = soilPhysical.KS[layer];

                            QTheta[layer, 1] -= flow / (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                            if (MathUtilities.IsGreaterThan(QTheta[layer, 1], soilPhysical.SAT[layer]))
                            {
                                flow = -(QTheta[layer, 1] - soilPhysical.SAT[layer]) * (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                                QTheta[layer, 1] = soilPhysical.SAT[layer];
                                Qpsi[layer, 1] = 0.0;
                                QK[layer, 1] = soilPhysical.KS[layer];

                                QTheta[layer, 0] -= flow / (QWeight[layer, 0] * soilPhysical.Thickness[layer]);
                                if (MathUtilities.IsGreaterThan(QTheta[layer, 0], soilPhysical.SAT[layer]))
                                {
                                    QTheta[layer, 0] = soilPhysical.SAT[layer];
                                    Qpsi[layer, 0] = 0.0;
                                    QK[layer, 0] = soilPhysical.KS[layer];
                                }
                                else
                                {
                                    Qpsi[layer, 0] = Suction(layer, QTheta[layer, 0]);
                                    QK[layer, 0] = SimpleK(layer, Qpsi[layer, 0]);
                                }
                            }
                            else
                            {
                                Qpsi[layer, 1] = Suction(layer, QTheta[layer, 1]);
                                QK[layer, 1] = SimpleK(layer, Qpsi[layer, 1]);
                            }
                        }
                        else
                        {
                            Qpsi[layer, 2] = Suction(layer, QTheta[layer, 2]);
                            QK[layer, 2] = SimpleK(layer, Qpsi[layer, 2]);
                        }
                    }
                    break;
                case 0:
                    if (flow > 0.0)
                    {
                        QTheta[layer, 1] += flow / (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsGreaterThan(QTheta[layer, 1], soilPhysical.SAT[layer]))
                        {
                            flow = 0.5 * (QTheta[layer, 1] - soilPhysical.SAT[layer]) * (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                            QTheta[layer, 1] = soilPhysical.SAT[layer];
                            Qpsi[layer, 1] = 0.0;
                            QK[layer, 1] = soilPhysical.KS[layer];

                            QTheta[layer, 0] += flow / (QWeight[layer, 0] * soilPhysical.Thickness[layer]);
                            if (MathUtilities.IsGreaterThan(QTheta[layer, 0], soilPhysical.SAT[layer]))
                            {
                                QTheta[layer, 0] = soilPhysical.SAT[layer];
                                Qpsi[layer, 0] = 0.0;
                                QK[layer, 0] = soilPhysical.KS[layer];
                            }
                            else
                            {
                                Qpsi[layer, 0] = Suction(layer, QTheta[layer, 0]);
                                QK[layer, 0] = SimpleK(layer, Qpsi[layer, 0]);
                            }

                            QTheta[layer, 2] += flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                            if (MathUtilities.IsGreaterThan(QTheta[layer, 2], soilPhysical.SAT[layer]))
                            {
                                QTheta[layer, 2] = soilPhysical.SAT[layer];
                                Qpsi[layer, 2] = 0.0;
                                QK[layer, 2] = soilPhysical.KS[layer];
                            }
                            else
                            {
                                Qpsi[layer, 2] = Suction(layer, QTheta[layer, 2]);
                                QK[layer, 2] = SimpleK(layer, Qpsi[layer, 2]);
                            }
                        }
                        else
                        {
                            Qpsi[layer, 1] = Suction(layer, QTheta[layer, 1]);
                            QK[layer, 1] = SimpleK(layer, Qpsi[layer, 1]);
                        }
                    }
                    else if (flow < 0.0)
                    {
                        // This flow shouldn't empty this section.
                        QTheta[layer, 1] += flow / (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsLessThan(QTheta[layer, 1], soilPhysical.AirDry[layer]))
                            System.Console.WriteLine("Too much water flow out from this section.");
                        Qpsi[layer, 1] = Suction(layer, QTheta[layer, 1]);
                        QK[layer, 1] = SimpleK(layer, Qpsi[layer, 1]);
                    }
                    break;
                case -1:
                    if (flow > 0.0)
                    {
                        QTheta[layer, 0] += flow / (QWeight[layer, 0] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsGreaterThan(QTheta[layer, 0], soilPhysical.SAT[layer]))
                        {
                            flow = (QTheta[layer, 0] - soilPhysical.SAT[layer]) * (QWeight[layer, 0] * soilPhysical.Thickness[layer]);
                            QTheta[layer, 0] = soilPhysical.SAT[layer];
                            Qpsi[layer, 0] = 0.0;
                            QK[layer, 0] = soilPhysical.KS[layer];

                            QTheta[layer, 1] += flow / (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                            if (MathUtilities.IsGreaterThan(QTheta[layer, 1], soilPhysical.SAT[layer]))
                            {
                                flow = (QTheta[layer, 1] - soilPhysical.SAT[layer]) * (QWeight[layer, 1] * soilPhysical.Thickness[layer]);
                                QTheta[layer, 1] = soilPhysical.SAT[layer];
                                Qpsi[layer, 1] = 0.0;
                                QK[layer, 1] = soilPhysical.KS[layer];

                                QTheta[layer, 2] += flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                                if (MathUtilities.IsGreaterThan(QTheta[layer, 2], soilPhysical.SAT[layer]))
                                {
                                    QTheta[layer, 2] = soilPhysical.SAT[layer];
                                    Qpsi[layer, 2] = 0.0;
                                    QK[layer, 2] = soilPhysical.KS[layer];
                                }
                                else
                                {
                                    Qpsi[layer, 2] = Suction(layer, QTheta[layer, 2]);
                                    QK[layer, 2] = SimpleK(layer, Qpsi[layer, 2]);
                                }
                            }
                            else
                            {
                                Qpsi[layer, 1] = Suction(layer, QTheta[layer, 1]);
                                QK[layer, 1] = SimpleK(layer, Qpsi[layer, 1]);
                            }
                        }
                        else
                        {
                            Qpsi[layer, 0] = Suction(layer, QTheta[layer, 0]);
                            QK[layer, 0] = SimpleK(layer, Qpsi[layer, 0]);
                        }
                    }
                    else if (flow < 0.0)
                    {
                        // This flow shouldn't empty this section.
                        QTheta[layer, 0] += flow / (QWeight[layer, 0] * soilPhysical.Thickness[layer]);
                        if (MathUtilities.IsLessThan(QTheta[layer, 0], soilPhysical.AirDry[layer]))
                            System.Console.WriteLine("Too much water flow out from this section.");
                        Qpsi[layer, 0] = Suction(layer, QTheta[layer, 0]);
                        QK[layer, 0] = SimpleK(layer, Qpsi[layer, 0]);
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
            // Current method allows the two sections to reach equilibrium easier than it should be.
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
                depth = (QDepth[layer + 1, 0] - QDepth[layer, 2]);
                psi_diff = (Qpsi[layer, 2] - Qpsi[layer + 1, 0]);
                grad = psi_diff / depth + 1.0;
                flow = (QK[layer, 2] + QK[layer + 1, 0]) / 2 * grad;
                if (Math.Abs(flow) < tolerance)
                    return flow;
                min_flow = 0.0;
                if (flow > 0.0)
                {
                    water_a = (QTheta[layer, 2] - SimpleTheta(layer, psiad)) * QWeight[layer, 2] * soilPhysical.Thickness[layer];
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
                    water_a = (QTheta[layer + 1, 0] - SimpleTheta(layer + 1, psiad)) * QWeight[layer + 1, 0] * soilPhysical.Thickness[layer + 1];
                    water_d = (soilPhysical.SAT[layer] - QTheta[layer, 2]) * QWeight[layer, 2] * soilPhysical.Thickness[layer];
                    max_flow = soilPhysical.KS[layer];
                    max_flow = -Math.Min(max_flow, Math.Min(water_a, water_d));
                    max_flow = Math.Max(max_flow, flow);
                    flow = max_flow;
                }

                old_flow = flow;
                for (int iteration = 0; iteration < 100; ++iteration)
                {
                    theta_1 = QTheta[layer, 2] - flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                    theta_2 = QTheta[layer + 1, 0] + flow / (QWeight[layer + 1, 0] * soilPhysical.Thickness[layer + 1]);
                    psi_1 = Suction(layer, theta_1);
                    psi_2 = Suction(layer + 1, theta_2);
                    grad = (psi_1 - psi_2) / depth + 1.0;

                    //if (MathUtilities.FloatsAreEqual(grad, 0.0))
                    //{
                    //    return flow;
                    //}

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

                theta = QTheta[layer, 2];
                flow = QK[layer, 2];
                water_a = (QTheta[layer, 2] - SimpleTheta(layer, psiad)) * QWeight[layer, 2] * soilPhysical.Thickness[layer];
                flow = Math.Min(flow, water_a);
                

                for (int i = 0; i < 100; ++i)
                {
                    if (flow <= flow_d)
                    {
                        flow_cum += flow;
                        return flow_cum;
                    }

                    flow = flow_d;
                    flow_cum += flow;
                    
                    theta -= flow / (QWeight[layer, 2] * soilPhysical.Thickness[layer]);
                    psi = Suction(layer, theta);
                    K = SimpleK(layer, psi);

                    flow = K;
                }

                return flow_cum;
            }
        }
    }
}
