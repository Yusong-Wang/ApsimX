﻿namespace Models.GrazPlan
{
    using APSIM.Shared.Utilities;
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Reflection;
    using System.Xml;

    /// <summary>
    /// Encapsulates a collection of stock genotype parameters. It can read the GrazPlan .prm
    /// files as well as the APSIM ruminant JSON file format.
    /// </summary>
    [Serializable]
    public class Genotypes
    {
        /// <summary>
        /// User supplied genotypes. These are searched first when looking for genotypes.
        /// </summary>
        private List<Genotype> genotypes = new List<Genotype>();

        /// <summary>Constructor.</summary>
        public Genotypes()
        {
            var resourceNames = Assembly.GetExecutingAssembly().GetManifestResourceNames().ToList()
                                        .Where(r => r.StartsWith("Models.Resources.GrazPlan.Genotypes."));
            foreach (var resourceName in resourceNames)
            {
                genotypes.Add(new Genotype(resourceName));
            }
        }

        /// <summary>Get a list of all genotypes.</summary>
        public IEnumerable<Genotype> All { get { return genotypes; } }

        /// <summary>Get a list of genotype names.</summary>
        public IEnumerable<string> Names { get { return All.Select(genotype => genotype.Name); } }

        /// <summary>
        /// Read a parameter set and append to the json array.
        /// </summary>
        /// <param name="xmlString">The XML string to read.</param>
        public void ReadPRM(string xmlString)
        {
            var xml = new XmlDocument();
            xml.LoadXml(xmlString);
            var parameters = xml.DocumentElement;

            ReadPRM(parameters);
        }

        /// <summary>Set the user specified genotypes.</summary>
        /// <param name="animalParameterSet">The user specified animal parameter set.</param>
        public void Add(AnimalParameterSet animalParameterSet)
        {
            Add(new Genotype(animalParameterSet));
        }

        /// <summary>Get a genotype. Throws if not found.</summary>
        /// <param name="genotypeName"></param>
        public Genotype Get(string genotypeName)
        {
            var foundGenotype = All.Where(genotype => genotype.Name.Equals(genotypeName, StringComparison.InvariantCultureIgnoreCase));
            if (foundGenotype.Count() == 0)
                throw new Exception($"Cannot find stock genotype {genotypeName}");
            return foundGenotype.First();
        }

        /// <summary>
        /// Read a parameter set and append to the json array.
        /// </summary>
        /// <param name="parameterNode">The XML parameter node to convert.</param>
        private void ReadPRM(XmlNode parameterNode)
        {
            Add(new Genotype(parameterNode));

            // recurse through child parameter sets.
            foreach (var child in XmlUtilities.ChildNodes(parameterNode, "set"))
                ReadPRM(child);
        }

        /// <summary>Add a genotype into the list of genotypes.</summary>
        /// <param name="genotypeToAdd">The genotype to add.</param>
        private void Add(Genotype genotypeToAdd)
        {
            var foundGenotype = genotypes.Find(genotype => genotype.Name == genotypeToAdd.Name);
            if (foundGenotype != null)
                genotypes.Remove(foundGenotype);
            genotypes.Add(genotypeToAdd);
        }
    }
}
