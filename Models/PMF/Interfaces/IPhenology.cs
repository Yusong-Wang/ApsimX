﻿using Models.Core;
using Models.PMF.Phen;

namespace Models.PMF.Interfaces
{
    /// <summary>
    /// An interface for a phenology model.
    /// </summary>
    /// <remarks>
    /// fixme - there's a lot of baggage here which should be removed.
    /// </remarks>
    public interface IPhenology : IModel
    {
        /// <summary>
        /// The current phenological phase.
        /// </summary>
        IPhase CurrentPhase { get; }

        /// <summary>
        /// A one based stage number.
        /// </summary>
        double Stage { get; set; }

        /// <summary>
        /// Gets the current zadok stage number. Used in manager scripts.
        /// </summary>
        double Zadok { get; }
    }
}
