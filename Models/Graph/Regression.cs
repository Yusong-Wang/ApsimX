namespace Models
{
    using APSIM.Shared.Utilities;
    using Models.Core;
    using Models.Storage;
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Drawing;
    using System.Globalization;
    using System.Linq;
    using System.Text;

    /// <summary>
    /// A regression model.
    /// </summary>
    [Serializable]
    [ViewName("UserInterface.Views.PropertyView")]
    [PresenterName("UserInterface.Presenters.PropertyPresenter")]
    [ValidParent(ParentType = typeof(Series))]
    [ValidParent(ParentType = typeof(Graph))]
    public class Regression : Model, ICachableGraphable
    {
        /// <summary>The stats from the regression</summary>
        private List<MathUtilities.RegrStats> stats = new List<MathUtilities.RegrStats>();

        /// <summary>The colours to use for each equation.</summary>
        private List<Color> equationColours = new List<Color>();

        /// <summary>
        /// Gets or sets a value indicating whether a regression should be shown for each series.
        /// </summary>
        /// <value><c>true</c> if [for each series]; otherwise, <c>false</c>.</value>
        [Description("Display regression line and equation for each series?")]
        public bool ForEachSeries { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether a regression should be shown for each series.
        /// </summary>
        /// <value><c>true</c> if [for each series]; otherwise, <c>false</c>.</value>
        [Description("Display 1:1 line?")]
        public bool showOneToOne { get; set; } = true;

        /// <summary>
        /// Gets or sets a value indicating whether a regression should be shown for each series.
        /// </summary>
        /// <value><c>true</c> if [for each series]; otherwise, <c>false</c>.</value>
        [Description("Display equation?")]
        public bool showEquation { get; set; } = true;

        /// <summary>Get a list of all actual series to put on the graph.</summary>
        /// <param name="storage">Storage service</param>
        /// <param name="simulationsFilter">Unused simulation names filter.</param>
        public IEnumerable<SeriesDefinition> GetSeriesDefinitions(IStorageReader storage, List<string> simulationsFilter = null)
        {
            Series seriesAncestor = FindAncestor<Series>();
            IEnumerable<SeriesDefinition> definitions;
            if (seriesAncestor == null)
            {
                Graph graph = FindAncestor<Graph>();
                if (graph == null)
                    throw new Exception("Regression model must be a descendant of a series");
                definitions = graph.FindAllChildren<Series>().SelectMany(s => s.GetSeriesDefinitions(storage, simulationsFilter));
            }
            else
                definitions = seriesAncestor.GetSeriesDefinitions(storage, simulationsFilter);

            return GetSeriesToPutOnGraph(storage, definitions, simulationsFilter);
        }

        /// <summary>Get a list of all actual series to put on the graph.</summary>
        /// <param name="storage">Storage service (required for access to checkpoint names).</param>
        /// <param name="definitions">Series definitions to be used (allows for caching of data).</param>
        /// <param name="simulationsFilter">Unused simulation names filter.</param>
        public IEnumerable<SeriesDefinition> GetSeriesToPutOnGraph(IStorageReader storage, IEnumerable<SeriesDefinition> definitions, List<string> simulationsFilter = null)
        {
            stats.Clear();
            equationColours.Clear();

            int checkpointNumber = 0;
            List<SeriesDefinition> regressionLines = new List<SeriesDefinition>();
            foreach (var checkpointName in storage.CheckpointNames)
            {
                // Get all x/y data
                List<double> x = new List<double>();
                List<double> y = new List<double>();
                foreach (SeriesDefinition definition in definitions)
                {
                    if (definition.CheckpointName == checkpointName)
                        if (definition.X is double[] && definition.Y is double[])
                        {
                            x.AddRange(definition.X as IEnumerable<double>);
                            y.AddRange(definition.Y as IEnumerable<double>);
                        }
                }

                if (ForEachSeries)
                {
                    // Display a regression line for each series.
                    // todo - should this also filter on checkpoint name?
                    int numDefinitions = definitions.Count();
                    foreach (SeriesDefinition definition in definitions)
                    {
                        if (definition.X is double[] && definition.Y is double[])
                        {
                            SeriesDefinition regressionSeries = PutRegressionLineOnGraph(definition.X, definition.Y, definition.Colour, null);
                            if (regressionSeries != null)
                            {
                                regressionLines.Add(regressionSeries);
                                equationColours.Add(definition.Colour);
                            }
                        }
                    }
                }
                else
                {
                    var regresionLineName = "Regression line";
                    if (checkpointName != "Current")
                        regresionLineName = "Regression line (" + checkpointName + ")";

                    // Display a single regression line for all data.
                    if (x.Count > 0 && y.Count == x.Count)
                    {
                        SeriesDefinition regressionSeries = PutRegressionLineOnGraph(x, y, ColourUtilities.ChooseColour(checkpointNumber), regresionLineName);
                        if (regressionSeries != null)
                        {
                            regressionLines.Add(regressionSeries);
                            equationColours.Add(ColourUtilities.ChooseColour(checkpointNumber));
                        }
                    }
                }

                if (showOneToOne)
                    regressionLines.Add(Put1To1LineOnGraph(x, y));

                checkpointNumber++;
            }

            return regressionLines;
        }
        
        /// <summary>Return a list of extra fields that the definition should read.</summary>
        /// <param name="seriesDefinition">The calling series definition.</param>
        /// <returns>A list of fields - never null.</returns>
        public IEnumerable<string> GetExtraFieldsToRead(SeriesDefinition seriesDefinition)
        {
            return new string[0];
        }

        /// <summary>Puts the regression line and 1:1 line on graph.</summary>
        /// <param name="x">The x data.</param>
        /// <param name="y">The y data.</param>
        /// <param name="colour">The colour of the regresion line.</param>
        /// <param name="title">The title to put in the legen.</param>
        private SeriesDefinition PutRegressionLineOnGraph(IEnumerable x, IEnumerable y, Color colour, string title)
        {
            MathUtilities.RegrStats stat = MathUtilities.CalcRegressionStats(title, y, x);
            if (stat != null)
            {
                stats.Add(stat);
                double minimumX = MathUtilities.Min(x);
                double maximumX = MathUtilities.Max(x);
                double minimumY = MathUtilities.Min(y);
                double maximumY = MathUtilities.Max(y);
                double lowestAxisScale = Math.Min(minimumX, minimumY);
                double largestAxisScale = Math.Max(maximumX, maximumY);

                var regressionDefinition = new SeriesDefinition
                    (title, colour,
                     new double[] { minimumX, maximumX },
                     new double[] { stat.Slope * minimumX + stat.Intercept, stat.Slope * maximumX + stat.Intercept });
                return regressionDefinition;
            }
            throw new Exception($"Unable to generate regression line for series {title} - there is no data");
        }

        /// <summary>Puts the 1:1 line on graph.</summary>
        /// <param name="x">The x data.</param>
        /// <param name="y">The y data.</param>
        private static SeriesDefinition Put1To1LineOnGraph(IEnumerable x, IEnumerable y)
        {
            IEnumerable<double> xValues = x.Cast<object>().Select(xi => xi is double ? (double)xi : Convert.ToDouble(xi, CultureInfo.InvariantCulture));
            IEnumerable<double> yValues = y.Cast<object>().Select(yi => yi is double ? (double)yi : Convert.ToDouble(yi, CultureInfo.InvariantCulture));
            MathUtilities.GetBounds(xValues, yValues, out double minX, out double maxX, out double minY, out double maxY);
            double lowestAxisScale = Math.Min(minX, minY);
            double largestAxisScale = Math.Max(maxX, maxY);

            return new SeriesDefinition
                ("1:1 line", Color.Empty,
                new double[] { lowestAxisScale, largestAxisScale },
                new double[] { lowestAxisScale, largestAxisScale },
                LineType.Dash, MarkerType.None);
        }

        /// <summary>Called by the graph presenter to get a list of all annotations to put on the graph.</summary>
        public IEnumerable<IAnnotation> GetAnnotations()
        {
            if (showEquation)
            {
                for (int i = 0; i < stats.Count; i++)
                {
                    // Add an equation annotation.
                    TextAnnotation equation = new TextAnnotation();
                    StringBuilder text = new StringBuilder();
                    text.AppendLine($"y = {stats[i].Slope:F2}x + {stats[i].Intercept:F2}, r2 = {stats[i].R2:F2}, n = {stats[i].n:F0}");
                    text.AppendLine($"NSE = {stats[i].NSE:F2}, ME = {stats[i].ME:F2}, MAE = {stats[i].MAE:F2}");
                    text.AppendLine($"RSR = {stats[i].RSR:F2}, RMSD = {stats[i].RMSE:F2}");
                    equation.Name = $"Regression{i}";
                    equation.text = text.ToString();
                    equation.colour = equationColours[i];
                    equation.leftAlign = true;
                    equation.textRotation = 0;
                    if (stats.Count > 1)
                    {
                        equation.x = double.MinValue;  // More than one stats equation. Use default positioning
                        equation.y = double.MinValue;
                    }
                    yield return equation;
                }
            }
        }
    }
}
