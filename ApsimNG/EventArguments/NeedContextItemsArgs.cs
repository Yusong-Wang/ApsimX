﻿namespace UserInterface.EventArguments
{
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Linq;
    using System.Reflection;
    using Models.Core;
    using System.Text;
    using APSIM.Shared.Utilities;
    using Intellisense;
    using System.Xml;
    using System.Drawing;
    using Models.Factorial;

    /// <summary>
    /// The editor view asks the presenter for context items. This structure
    /// is used to do that
    /// </summary>
    public class NeedContextItemsArgs : EventArgs
    {
        /// <summary>
        /// The name of the object that needs context items.
        /// </summary>
        public string ObjectName;

        /// <summary>
        /// The items returned from the presenter back to the view
        /// </summary>
        public List<ContextItem> AllItems;

        /// <summary>
        /// Context item information
        /// </summary>
        public List<string> Items;
#if NETFRAMEWORK
        /// <summary>
        /// Completion data.
        /// </summary>
        public List<CompletionData> CompletionData { get; set; }
#endif
        /// <summary>
        /// Co-ordinates at which the intellisense window should be displayed.
        /// </summary>
        public Point Coordinates { get; set; }

        /// <summary>
        /// Source code for which we need completion options.
        /// </summary>
        public string Code { get; set; }

        /// <summary>
        /// Offset of the caret in the source code.
        /// </summary>
        public int Offset { get; set; }

        /// <summary>
        /// True iff this intellisense request was generated by the user pressing control space.
        /// </summary>
        public bool ControlSpace { get; set; }

        /// <summary>
        /// True iff this intellisense request was generated by the user pressing control space.
        /// </summary>
        public bool ControlShiftSpace { get; set; }

        /// <summary>
        /// The line that the caret is on.
        /// </summary>
        public int LineNo { get; set; }

        /// <summary>
        /// The column that the caret is on.
        /// </summary>
        public int ColNo { get; set; }

        /// <summary>
        /// A dictionary mapping assemblies to their xml documentation.
        /// Used to cache the xml documents, speeding up intellisense operations.
        /// </summary>
        private static Dictionary<string, XmlDocument> documentation = new Dictionary<string, XmlDocument>();

        /// <summary>
        /// The view is asking for variable names for its intellisense.
        /// </summary>
        /// <param name="atype">Data type for which we want completion options.</param>
        /// <param name="properties">If true, property suggestions will be generated.</param>
        /// <param name="methods">If true, method suggestions will be generated.</param>
        /// <param name="publishedEvents">If true, published events will be returned.</param>
        /// <param name="subscribedEvents">If true, subscribed events will be returned.</param>
        /// <returns>List of completion options.</returns>
        public static List<ContextItem> ExamineTypeForContextItems(Type atype, bool properties, bool methods, bool publishedEvents, bool subscribedEvents)
        {
            List<ContextItem> allItems = new List<ContextItem>();

            // find the properties and methods
            if (atype != null)
            {
                if (properties)
                {
                    foreach (PropertyInfo property in atype.GetProperties(BindingFlags.Instance | BindingFlags.Public))
                    {
                        VariableProperty var = new VariableProperty(atype, property);
                        ContextItem item = new ContextItem
                        {
                            Name = var.Name,
                            IsProperty = true,
                            IsEvent = false,
                            IsWriteable = !var.IsReadOnly,
                            TypeName = GetTypeName(var.DataType),
                            Descr = GetDescription(property),
                            Units = var.Units
                        };
                        allItems.Add(item);
                    }
                }

                if (methods)
                {
                    foreach (MethodInfo method in atype.GetMethods(BindingFlags.Instance | BindingFlags.Public | BindingFlags.FlattenHierarchy))
                    {
                        if (!method.Name.StartsWith("get_") && !method.Name.StartsWith("set_") &&
                            !method.Name.StartsWith("add_") && !method.Name.StartsWith("remove_"))
                        {
                            ContextItem item = new ContextItem
                            {
                                Name = method.Name,
                                IsProperty = false,
                                IsEvent = false,
                                IsMethod = true,
                                IsWriteable = false,
                                TypeName = method.ReturnType.Name,
                                Descr = GetDescription(method),
                                Units = string.Empty
                            };

                            // build a parameter string representation
                            ParameterInfo[] allparams = method.GetParameters();
                            StringBuilder paramText = new StringBuilder("( ");
                            if (allparams.Count() > 0)
                            {
                                for (int p = 0; p < allparams.Count(); p++)
                                {
                                    ParameterInfo parameter = allparams[p];
                                    paramText.Append(parameter.ParameterType.Name + " " + parameter.Name);
                                    if (parameter.DefaultValue != DBNull.Value)
                                        paramText.Append(" = " + parameter.DefaultValue?.ToString() ?? "null");
                                    if (p < allparams.Count() - 1)
                                        paramText.Append(", ");
                                }
                            }
                            paramText.Append(" )");
                            item.ParamString = paramText.ToString();

                            allItems.Add(item);
                        }
                    }
                }

                if (publishedEvents)
                {
                    foreach (EventInfo evnt in atype.GetEvents(BindingFlags.Instance | BindingFlags.Public))
                    {
                        NeedContextItemsArgs.ContextItem item = new NeedContextItemsArgs.ContextItem();
                        item.Name = evnt.Name;
                        item.IsProperty = true;
                        item.IsEvent = true;
                        item.IsMethod = false;
                        item.IsWriteable = false;
                        item.TypeName = evnt.ReflectedType.Name;
                        item.Descr = GetDescription(evnt);
                        item.Units = "";
                        allItems.Add(item);
                    }
                }

                if (subscribedEvents)
                {
                    foreach (MethodInfo method in atype.GetMethods(/*BindingFlags.Instance |*/ BindingFlags.NonPublic | BindingFlags.FlattenHierarchy))
                    {
                        EventSubscribeAttribute subscriberAttribute = (EventSubscribeAttribute)ReflectionUtilities.GetAttribute(method, typeof(EventSubscribeAttribute), false);
                        if (subscriberAttribute != null)
                        {
                            NeedContextItemsArgs.ContextItem item = new NeedContextItemsArgs.ContextItem();
                            item.Name = subscriberAttribute.ToString();
                            item.IsProperty = false;
                            item.IsEvent = false;
                            item.IsMethod = true;
                            item.IsWriteable = false;
                            item.TypeName = atype.Name;
                            item.Descr = GetDescription(atype);
                            item.Units = "";
                            allItems.Add(item);
                        }
                    }
                }
            }

            allItems.Sort(delegate(ContextItem c1, ContextItem c2) { return c1.Name.CompareTo(c2.Name); });
            return allItems;
        }

        /// <summary>
        /// The view is asking for variable names for its intellisense.
        /// </summary>
        /// <param name="o">Fully- or partially-qualified object name for which we want completion options.</param>
        /// <param name="properties">If true, property suggestions will be generated.</param>
        /// <param name="methods">If true, method suggestions will be generated.</param>
        /// <param name="publishedEvents">If true, published events will be returned.</param>
        /// <param name="subscribedEvents">If true, subscribed events will be returned.</param>
        /// <returns>List of completion options.</returns>
        private static List<ContextItem> ExamineObjectForContextItems(object o, bool properties, bool methods, bool publishedEvents, bool subscribedEvents)
        {
            List<ContextItem> allItems;
            Type objectType = o is Type ? o as Type : o.GetType();
            allItems = ExamineTypeForContextItems(objectType, properties, methods, publishedEvents, subscribedEvents);
            
            // add in the child models.
            if (o is IModel)
            {
                foreach (IModel model in (o as IModel).Children)
                {
                    if (allItems.Find(m => m.Name == model.Name) == null)
                    {
                        NeedContextItemsArgs.ContextItem item = new NeedContextItemsArgs.ContextItem();
                        item.Name = model.Name;
                        item.IsProperty = false;
                        item.IsEvent = false;
                        item.IsWriteable = false;
                        item.TypeName = model.GetType().Name;
                        item.Units = string.Empty;
                        allItems.Add(item);
                    }
                }
                allItems.Sort(delegate(ContextItem c1, ContextItem c2) { return c1.Name.CompareTo(c2.Name); });
            }
            return allItems;
        }

        /// <summary>
        /// The view is asking for variable names.
        /// </summary>
        /// <param name="relativeTo">Model in the simulations tree which owns the editor.</param>
        /// <param name="objectName">Fully- or partially-qualified object name for which we want completion options.</param>
        /// <param name="properties">If true, property suggestions will be generated.</param>
        /// <param name="methods">If true, method suggestions will be generated.</param>
        /// <param name="events">If true, event suggestions will be generated.</param>
        /// <returns>List of completion options.</returns>
        public static List<ContextItem> ExamineModelForNames(IModel relativeTo, string objectName, bool properties, bool methods, bool events)
        {
            // TODO : refactor cultivar and report activity ledger presenters so they use the intellisense presenter. 
            // These are the only two presenters which still use this intellisense method.
            if (objectName == string.Empty)
                objectName = ".";

            object o = null;
            IModel replacementModel = relativeTo.FindByPath(".Simulations.Replacements")?.Value as IModel;
            if (replacementModel != null)
            {
                try
                {
                    o = replacementModel.FindByPath(objectName)?.Value as IModel;
                }
                catch (Exception) {  }
            }

            if (o == null)
            {
                try
                {
                    o = relativeTo.FindByPath(objectName)?.Value;
                }
                catch (Exception) { }
            }
            
            if (o == null && relativeTo.Parent is Replacements)
            {
                // Model 'relativeTo' could be under replacements. Look for the first simulation and try that.
                IModel simulation = relativeTo.Parent.Parent.FindInScope<Simulation>();
                try
                {
                    o = simulation.FindByPath(objectName)?.Value as IModel;
                }
                catch (Exception) { }
            }

            if (o != null)
            {
                return ExamineObjectForContextItems(o, properties, methods, events, false);
            }

            return new List<ContextItem>();
        }

        /// <summary>
        /// Generates a list of context items for given model.
        /// Uses <see cref="GetNodeFromPath(Model, string)"/> to get the model reference.
        /// </summary>
        /// <param name="relativeTo">Model that the string is relative to.</param>
        /// <param name="objectName">Name of the model that we want context items for.</param>
        /// <param name="properties">Search for properties of the model?</param>
        /// <param name="methods">Search for methods of the model?</param>
        /// <param name="publishedEvents">If true, published events will be returned.</param>
        /// <param name="subscribedEvents">If true, subscribed events will be returned.</param>
        public static List<ContextItem> ExamineModelForContextItemsV2(Model relativeTo, string objectName, bool properties, bool methods, bool publishedEvents, bool subscribedEvents)
        {
            List<ContextItem> contextItems = new List<ContextItem>();
            object node = GetNodeFromPath(relativeTo, objectName);
            if (node == null)
                node = relativeTo.FindByPath(objectName)?.Value;
            if (node != null)
            {
                contextItems = ExamineObjectForContextItems(node, properties, methods, publishedEvents, subscribedEvents);
            }
            return contextItems;
        }

        public static MethodInfo GetMethodInfo(Model relativeTo, string methodName, string objectName)
        {
            object node = GetNodeFromPath(relativeTo, objectName);
            if (node != null)
                return node.GetType().GetMethod(methodName);
            return null;
        }

        /// <summary>
        /// A new method for finding a model/object from a path in the simulations tree.
        /// Finds the node (whose name is surrounded by square brackets). From there, it looks for each
        /// successive period-delimited child or property given in the path string.
        /// </summary>
        /// <param name="relativeTo">Object in the simulations tree.</param>
        /// <param name="objectName">Name of the object or model for which we want completion options.</param>
        /// <returns></returns>
        private static object GetNodeFromPath(Model relativeTo, string objectName)
        {       
            string modelNamePattern = @"\[[A-Za-z\s]+[A-Za-z0-9\s_]*\]";
            object node = null;
            var matches = System.Text.RegularExpressions.Regex.Matches(objectName, modelNamePattern);
            if (matches.Count <= 0)
            {
                // object name doesn't contain square brackets.
                string textBeforeFirstDot = objectName;
                if (objectName.Contains("."))
                    if (objectName.StartsWith("."))
                        textBeforeFirstDot = textBeforeFirstDot.Substring(1, textBeforeFirstDot.Length - 1);
                    else
                        textBeforeFirstDot = textBeforeFirstDot.Substring(0, textBeforeFirstDot.IndexOf('.'));
                node = relativeTo.FindInScope(textBeforeFirstDot);
                if (node == null)
                    node = relativeTo.FindByPath(objectName)?.Value;
            }
            else
            {
                // Get the raw model name without square brackets.
                string modelName = matches[0].Value.Replace("[", "").Replace("]", "");

                // Get the node in the simulations tree corresponding to the model name which was surrounded by square brackets.
                node = relativeTo.FindInScope(modelName);

                // If we're under replacements we won't be able to find some simulation-
                // related nodes such as weather/soil/etc. In this scenario, we should
                // search through all models, not just those in scope.
                if (node == null && relativeTo.FindAncestor<Replacements>() != null)
                {
                    node = relativeTo.FindAncestor<Simulations>().FindAllDescendants().FirstOrDefault(child => child.Name == modelName);

                    // If we still failed, try a lookup on type name.
                    if (node == null)
                        node = relativeTo.FindAncestor<Simulations>().FindAllDescendants().FirstOrDefault(x => x.GetType().Name == modelName);
                }

                if (node == null && relativeTo.FindAncestor<Factors>() != null)
                {
                    relativeTo = relativeTo.FindAncestor<Experiment>();
                    if (relativeTo != null)
                    {
                        node = relativeTo.FindInScope(modelName);
                        if (node == null)
                            node = relativeTo.FindAllDescendants().FirstOrDefault(x => x.GetType().Name == modelName);
                    }
                    if (node == null)
                        return null;
                }
            }

            // If the object name string does not contain any children/properties 
            // (e.g. it doesn't contain any periods), we can return immediately.
            if (!objectName.Contains("."))
                return node;

            objectName = objectName.Substring(objectName.IndexOf('.') + 1);

            // Iterate over the 'child' models/properties.
            // childName is the next child we're looking for. e.g. in "[Wheat].Leaf", the first childName will be "Leaf".
            string[] namePathBits = StringUtilities.SplitStringHonouringBrackets(objectName, ".", '[', ']');
            for (int i = 0; i < namePathBits.Length; i++)
            {
                if (node == null)
                    return null;
                string childName = namePathBits[i];

                int squareBracketIndex = childName.IndexOf('[');
                if (squareBracketIndex == 0)
                {
                    // User has typed something like [Wheat].[...]
                    throw new Exception("Unable to parse child or property " + childName);
                }
                if (squareBracketIndex > 0) // childName contains square brackets - it may be an IList element
                    childName = childName.Substring(0, squareBracketIndex);
                
                // First, check the child models. 
                if (node is IModel)
                    node = (node as IModel).Children.FirstOrDefault(c => c.Name == childName) ?? node;

                // If we couldn't find a matching child model, we check the model/object's properties.

                // This expression evaluates to true if node is not an IModel.
                if ((node as IModel)?.Name != childName)
                {
                    // Node cannot be null here.
                    try
                    {
                        Type propertyType = node is Type ? node as Type : node.GetType();
                        PropertyInfo property = propertyType.GetProperties(BindingFlags.Public | BindingFlags.Instance).FirstOrDefault(p => p.Name == childName);

                        // If we couldn't find any matching child models or properties, all we can do is return.
                        if (property == null)
                            return null;

                        // Try to set node to the value of the property.
                        try
                        {
                            node = ReflectionUtilities.GetValueOfFieldOrProperty(childName, node);
                        }
                        catch (TargetInvocationException err)
                        {
                            if (err.InnerException is NullReferenceException)
                                // Some properties depend on links being resolved which is not
                                // always the case (ie in an intellisense context).
                                node = null;
                            else
                                throw;
                        }

                        if (node == null)
                        {
                            // This property has the correct name. If the property's type provides a parameterless constructor, we can use 
                            // reflection to instantiate an object of that type and assign it to the node variable. 
                            // Otherwise, we assign the type itself to node.
                            if (property.PropertyType.GetConstructor(Type.EmptyTypes) == null)
                                node = property.PropertyType;
                            else
                                node = Activator.CreateInstance(property.PropertyType);
                        }
                    }
                    catch
                    {
                        // Any 'real' errors should be displayed to the user, so they should be caught 
                        // in a presenter which can access the explorer presenter.
                        // Because of this, any unhandled exceptions here will kill the intellisense 
                        // generation operation, and we still have a few tricks up our sleeve.
                        return null;
                    }
                }

                if (squareBracketIndex > 0)
                {
                    // We have found the node, but the node is an IList of some sort, and we are actually interested in a specific element.

                    int closingBracketIndex = namePathBits[i].IndexOf(']');
                    if (closingBracketIndex <= 0 || (closingBracketIndex - squareBracketIndex) < 1)
                        return null;

                    string textBetweenBrackets = namePathBits[i].Substring(squareBracketIndex + 1, closingBracketIndex - squareBracketIndex - 1);
                    if (node is IList)
                    {
                        int index = -1;
                        if (Int32.TryParse(textBetweenBrackets, out index))
                        {
                            IList nodeList = node as IList;
                            if (index > nodeList.Count || index <= 0)
                                node = node.GetType().GetInterfaces().Where(x => x.IsGenericType && x.GetGenericTypeDefinition() == typeof(IEnumerable<>)).Select(x => x.GetGenericArguments()[0]).FirstOrDefault();
                            else
                                node = nodeList[index - 1];
                        }
                        else
                            throw new Exception("Unable to access element \"" + textBetweenBrackets + "\" of list \"" + namePathBits[i] + "\"");
                    }
                    else if (node is IDictionary)
                    {
                        node = (node as IDictionary)[textBetweenBrackets];
                    }
                    else
                        throw new Exception("Unable to parse child or property name " + namePathBits[i]);
                }
                squareBracketIndex = -1;
            }
             return node;
        }

        private static string GetTypeName(Type type)
        {
            if (type.IsGenericType)
            {
                string name = type.Name;
                if (name.Contains("`"))
                    name = name.Substring(0, name.IndexOf('`'));
                return $"{name}<{string.Join(", ", type.GenericTypeArguments.Select(t => GetTypeName(t)))}>";
            }
            return type.Name;
        }

        /// <summary>
        /// Gets the contents of a property's summary tag, or, if the summary tag doesn't exist,
        /// a <see cref="DescriptionAttribute"/>.
        /// </summary>
        /// <param name="member">Property whose documentation will be retrieved.</param>
        /// <returns>
        /// Contents of a summary tag (if available), or a description attribute,
        /// or an empty string if neither of these are available.
        /// </returns>
        public static string GetDescription(MemberInfo member)
        {
            if (member == null)
                return string.Empty;

            // The member's documentation doesn't reside in the compiled assembly - it's actually stored in
            // an xml documentation file which usually sits next to the assembly on disk.
            string documentationFile = System.IO.Path.ChangeExtension(member.Module.Assembly.Location, "xml");
            if (!System.IO.File.Exists(documentationFile))
            {
                // If the documentation file doesn't exist, this member is probably a member of a manager script.
                // These members usually have a description attribute which can be used instead.
                Attribute description = member.GetCustomAttribute(typeof(DescriptionAttribute));
                // If the property has no description attribute, just return an empty string.
                return description == null ? string.Empty : description.ToString();
            }

            XmlDocument doc = new XmlDocument();
            if (documentation.ContainsKey(documentationFile))
                doc = documentation[documentationFile];
            else
            {
                doc.Load(documentationFile);
                documentation.Add(documentationFile, doc);
            }
            string memberPrefix = string.Empty;
            if (member is PropertyInfo)
                memberPrefix = "P";
            else if (member is FieldInfo)
                // This shouldn't be run for public fields, but it doesn't hurt to be thorough.
                memberPrefix = "F";
            else if (member is MethodInfo)
                memberPrefix = "M";
            else if (member is EventInfo)
                memberPrefix = "E";

            string xPath = string.Format("//member[@name='{0}:{1}.{2}']/summary[1]", memberPrefix, member.DeclaringType.FullName, member.Name);
            XmlNode summaryNode = doc.SelectSingleNode(xPath);
            if (summaryNode == null || summaryNode.ChildNodes.Count < 1)
                return string.Empty;
            // The summary tag often spans multiple lines, which means that the text inside usually
            // starts and ends with newlines (and indentation), which we don't want to display.
            return summaryNode.InnerText.Trim(Environment.NewLine.ToCharArray()).Trim();
        }

        /// <summary>
        /// Gets the XML documentation for a particular parameter of a method.
        /// </summary>
        /// <param name="method">The method.</param>
        /// <param name="parameterName">The parameter.</param>
        public static string GetDescription(MethodInfo method, string parameterName)
        {
            if (method == null)
                return string.Empty;

            // The method's documentation doesn't reside in the compiled assembly - it's actually stored in
            // an xml documentation file which usually sits next to the assembly on disk.
            string documentationFile = System.IO.Path.ChangeExtension(method.Module.Assembly.Location, "xml");
            if (!System.IO.File.Exists(documentationFile))
                return string.Empty;

            XmlDocument doc = new XmlDocument();

            if (documentation.ContainsKey(documentationFile))
                doc = documentation[documentationFile];
            else
            {
                doc.Load(documentationFile);
                documentation.Add(documentationFile, doc);
            }

            string xPath = string.Format("//member[starts-with(@name, 'M:{0}.{1}')]/param[@name='{2}']", method.DeclaringType.FullName, method.Name, parameterName);
            XmlNode paramNode = doc.SelectSingleNode(xPath);
            if (paramNode == null || paramNode.ChildNodes.Count < 1)
                return string.Empty;
            // The summary tag often spans multiple lines, which means that the text inside usually
            // starts and ends with newlines (and indentation), which we don't want to display.
            return paramNode.InnerText.Trim(Environment.NewLine.ToCharArray()).Trim();
        }

        /// <summary>
        /// Complete context item information
        /// </summary>
        public class ContextItem
        {
            /// <summary>
            /// Name of the item
            /// </summary>
            public string Name;

            /// <summary>
            /// The return type as a string
            /// </summary>
            public string TypeName;

            /// <summary>
            /// Units string
            /// </summary>
            public string Units;

            /// <summary>
            /// The description string
            /// </summary>
            public string Descr;

            /// <summary>
            /// This is an event.
            /// </summary>
            public bool IsEvent;

            /// <summary>
            /// This is a method.
            /// </summary>
            public bool IsMethod;

            /// <summary>
            /// String that represents the parameter list
            /// </summary>
            public string ParamString;

            /// <summary>
            /// This is a property
            /// </summary>
            public bool IsProperty;

            /// <summary>
            /// This property is writeable
            /// </summary>
            public bool IsWriteable;

            /// <summary>
            /// The property is a child model.
            /// </summary>
            public bool IsChildModel;

            /// <summary>
            /// Returns a new instance of a context item representing
            /// an event.
            /// </summary>
            /// <param name="name">Name of the event.</param>
            /// <param name="description">Description of the event.</param>
            /// <param name="eventType">Event handler type.</param>
            public static ContextItem NewEvent(string name, string description, Type eventType)
            {
                return new ContextItem()
                {
                    IsChildModel = false,
                    IsEvent = true,
                    IsMethod = false,
                    IsProperty = false,
                    IsWriteable = false,
                    Name = name,
                    Descr = description,
                    TypeName = eventType.Name,
                };
            }
        }
    } 
}
