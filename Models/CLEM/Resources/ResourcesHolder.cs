﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;  //enumerator
using Newtonsoft.Json;
using System.Runtime.Serialization;
using Models.Core;
using Models.CLEM.Activities;
using Models.CLEM.Groupings;
using System.ComponentModel.DataAnnotations;
using Models.Core.Attributes;
using APSIM.Shared.Utilities;

namespace Models.CLEM.Resources
{
    ///<summary>
    /// Manger for all resources available to the model
    ///</summary> 
    [Serializable]
    [ViewName("UserInterface.Views.GridView")]
    [PresenterName("UserInterface.Presenters.PropertyPresenter")]
    [ValidParent(ParentType = typeof(ZoneCLEM))]
    [ValidParent(ParentType = typeof(Market))]
    [Description("This holds all resource groups used in the CLEM simulation")]
    [Version(1, 0, 1, "")]
    [HelpUri(@"Content/Features/Resources/ResourcesHolder.htm")]
    public class ResourcesHolder: CLEMModel, IValidatableObject
    {
        /// <summary>
        /// List of the all the Resource Groups.
        /// </summary>
        [JsonIgnore]
        private List<IModel> ResourceGroupList;

        private void InitialiseResourceGroupList()
        {
            if(ResourceGroupList == null)
            {
                ResourceGroupList = this.FindAllChildren<IModel>().Where(a => a.Enabled).ToList();
            }
        }

        private IModel GetGroupByName(string name)
        {
            InitialiseResourceGroupList();
            return ResourceGroupList.Find(x => x.Name == name);
        }

        private IModel GetGroupByType(Type type)
        {
            InitialiseResourceGroupList();
            return ResourceGroupList.Find(x => x.GetType() == type);
        }

        /// <summary>
        /// Finds a shared marketplace
        /// </summary>
        /// <returns>Market</returns>
        [JsonIgnore]
        public Market FoundMarket { get; private set; }

        /// <summary>
        /// Determines if a market has been located
        /// </summary>
        /// <returns>True or false</returns>
        public bool MarketPresent { get { return !(FoundMarket is null); } }

        /// <summary>
        /// Determines whether resource items of the specified group type exist 
        /// </summary>
        /// <param name="resourceGroupType"></param>
        /// <returns></returns>
        public bool ResourceItemsExist(Type resourceGroupType)
        {
            Model resourceGroup = this.FindAllChildren().Where(c => resourceGroupType.IsAssignableFrom(c.GetType())).FirstOrDefault() as Model;
            if (resourceGroup != null && resourceGroup.Enabled)
            {
                return resourceGroup.Children.Where(a => a.Enabled).Count() > 0;
            }
            return false;
        }

        /// <summary>
        /// Determines whether resource group of the specified type exist 
        /// </summary>
        /// <param name="resourceGroupType"></param>
        /// <returns></returns>
        public bool ResourceGroupExist(Type resourceGroupType)
        {
            Model resourceGroup = this.FindAllChildren().Where(c => resourceGroupType.IsAssignableFrom(c.GetType())).FirstOrDefault() as Model;
            return (resourceGroup != null && resourceGroup.Enabled);
        }

        /// <summary>
        /// Get resource by name
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public object GetResourceGroupByName(string name)
        {
            InitialiseResourceGroupList();
            return ResourceGroupList.Find(x => x.Name == name);
        }

        /// <summary>
        /// Get resource by type
        /// </summary>
        /// <param name="resourceGroupType"></param>
        /// <returns></returns>
        public object GetResourceGroupByType(Type resourceGroupType)
        {
            InitialiseResourceGroupList();
            return ResourceGroupList.Find(x => x.GetType() == resourceGroupType);
        }

        /// <summary>
        /// Retrieve a ResourceType from a ResourceGroup based on a request item including filter and sort options
        /// </summary>
        /// <param name="request">A resource request item</param>
        /// <param name="missingResourceAction">Action to take if requested resource group not found</param>
        /// <param name="missingResourceTypeAction">Action to take if requested resource type not found</param>
        /// <returns>A reference to the item of type Model</returns>
        public IModel GetResourceItem(ResourceRequest request, OnMissingResourceActionTypes missingResourceAction, OnMissingResourceActionTypes missingResourceTypeAction)
        {
            if (request.FilterDetails != null)
            {
                if (request.ResourceType == null)
                {
                    string errorMsg = String.Format("Resource type must be supplied in resource request from [a={0}]", request.ActivityModel.Name);
                    throw new Exception(errorMsg);
                }

                IModel resourceGroup = this.GetGroupByType(request.ResourceType);
                if(resourceGroup== null)
                {
                    string errorMsg = String.Format("@error:Unable to locate resources of type [r{0}] for [a={1}]", request.ResourceType, request.ActivityModel.Name);
                    switch (missingResourceAction)
                    {
                        case OnMissingResourceActionTypes.ReportErrorAndStop:
                            throw new Exception(errorMsg);
                        case OnMissingResourceActionTypes.ReportWarning:
                            errorMsg = errorMsg.Replace("@error:", "");
                            Summary.WriteWarning(request.ActivityModel, errorMsg);
                            break;
                        default:
                            break;
                    }
                    return null;
                }

                // get list of children matching the conditions in filter
                // and return the lowest item that has enough time available
                object resourceGroupObject = resourceGroup as object;
                switch (resourceGroupObject.GetType().ToString())
                {
                    case "Models.CLEM.Resources.Labour":
                        // get matching labour types
                        // use activity uid to ensure unique for this request
                        List<LabourType> items = (resourceGroup as Labour).Items;
                        items = items.Filter(request.FilterDetails.FirstOrDefault() as Model);
                        items = items.Where(a => a.LastActivityRequestID != request.ActivityID).ToList();
                        if (items.Where(a => a.Amount >= request.Required).Count()>0)
                        {
                            // get labour least available but with the amount needed
                            return items.Where(a => a.Amount >= request.Required).OrderByDescending(a => a.Amount).FirstOrDefault();
                        }
                        else
                        {
                            // get labour with most available but with less than the amount needed
                            return items.OrderByDescending(a => a.Amount).FirstOrDefault();
                        }
                    default:
                        string errorMsg = "Resource cannot be filtered. Filtering not implemented for [r=" + resourceGroupObject.GetType().ToString() + "] from activity [a=" + request.ActivityModel.Name + "]";
                        Summary.WriteWarning(request.ActivityModel, errorMsg);
                        throw new Exception(errorMsg);
                }
            }
            else
            {
                // check style of ResourceTypeName used
                // this is either "Group.Type" from dropdown menus or "Type" only. 
                if (request.ResourceTypeName.Contains("."))
                {
                    return GetResourceItem(request.ActivityModel, request.ResourceTypeName, missingResourceAction, missingResourceTypeAction);
                }
                else
                {
                    return GetResourceItem(request.ActivityModel, request.ResourceType, request.ResourceTypeName, missingResourceAction, missingResourceTypeAction);
                }
            }
        }

        /// <summary>
        /// Retrieve a ResourceType from a ResourceGroup with specified names
        /// </summary>
        /// <param name="requestingModel">name of model requesting resource</param>
        /// <param name="resourceGroupType">Type of the resource group</param>
        /// <param name="resourceItemName">Name of the resource item</param>
        /// <param name="missingResourceAction">Action to take if requested resource group not found</param>
        /// <param name="missingResourceTypeAction">Action to take if requested resource type not found</param>
        /// <returns>A reference to the item of type object</returns>
        public IModel GetResourceItem(Model requestingModel, Type resourceGroupType, string resourceItemName, OnMissingResourceActionTypes missingResourceAction, OnMissingResourceActionTypes missingResourceTypeAction)
        {
            // locate specified resource
            Model resourceGroup = this.FindAllChildren().Where(c => resourceGroupType.IsAssignableFrom(c.GetType())).FirstOrDefault() as Model;
            if (resourceGroup != null)
            {
                IModel resource = resourceGroup.Children.Where(a => a.Name == resourceItemName & a.Enabled).FirstOrDefault();
                if (resource == null)
                {
                    string errorMsg = String.Format("@error:Unable to locate resources item [r={0}] in resources [r={1}] for [a={2}]", resourceItemName, resourceGroupType.ToString(), requestingModel.Name);
                    switch (missingResourceTypeAction)
                    {
                        case OnMissingResourceActionTypes.ReportErrorAndStop:
                            throw new Exception(errorMsg);
                        case OnMissingResourceActionTypes.ReportWarning:
                            Summary.WriteWarning(requestingModel, errorMsg);
                            break;
                        default:
                            break;
                    }
                    return null;
                }
                return resource;
            }
            else
            {
                string errorMsg = String.Format("@error:Unable to locate resources of type [r={0}] for [a={1}]", resourceGroupType.ToString(), requestingModel.Name);
                switch (missingResourceAction)
                {
                    case OnMissingResourceActionTypes.ReportErrorAndStop:
                        throw new Exception(errorMsg);
                    case OnMissingResourceActionTypes.ReportWarning:
                        errorMsg = errorMsg.Replace("@error:", "");
                        Summary.WriteWarning(requestingModel, errorMsg);
                        break;
                    default:
                        break;
                }
                return null;
            }
        }

        /// <summary>
        /// Retrieve a ResourceType from a ResourceGroup with specified names
        /// </summary>
        /// <param name="requestingModel">name of model requesting resource</param>
        /// <param name="resourceGroupAndItem">Period separated list of resource group and type</param>
        /// <param name="missingResourceAction">Action to take if requested resource group not found</param>
        /// <param name="missingResourceTypeAction">Action to take if requested resource type not found</param>
        /// <returns>A reference to the item of type object</returns>
        public IModel GetResourceItem(Model requestingModel, string resourceGroupAndItem, OnMissingResourceActionTypes missingResourceAction, OnMissingResourceActionTypes missingResourceTypeAction)
        {
            if(resourceGroupAndItem == null)
            {
                resourceGroupAndItem = " . ";
            }

            // locate specified resource
            string[] names = resourceGroupAndItem.Split('.');
            if(names.Count()!=2)
            {
                string errorMsg = String.Format("@error:Invalid resource group and type string for [{0}], expecting 'ResourceName.ResourceTypeName'. Value provided [{1}] ", requestingModel.Name, resourceGroupAndItem);
                throw new Exception(errorMsg);
            }

            Model resourceGroup = this.GetResourceGroupByName(names[0]) as Model;
            if (resourceGroup != null)
            {
                IModel resource = resourceGroup.Children.Where(a => a.Name == names[1] & a.Enabled).FirstOrDefault();
                if (resource == null)
                {
                    string errorMsg = String.Format("@error:Unable to locate resources item [r={0}] in resources [r={1}] for [a={2}]", names[1], names[0], requestingModel.Name);
                    switch (missingResourceTypeAction)
                    {
                        case OnMissingResourceActionTypes.ReportErrorAndStop:
                            throw new Exception(errorMsg);
                        case OnMissingResourceActionTypes.ReportWarning:
                            Summary.WriteWarning(requestingModel, errorMsg);
                            break;
                        default:
                            break;
                    }
                    return null;
                }
                return resource;
            }
            else
            {
                string errorMsg = String.Format("@error:Unable to locate resources of type [r={0}] for [a={1}]", names[0], requestingModel.Name);
                switch (missingResourceAction)
                {
                    case OnMissingResourceActionTypes.ReportErrorAndStop:
                        throw new Exception(errorMsg);
                    case OnMissingResourceActionTypes.ReportWarning:
                        errorMsg = errorMsg.Replace("@error:", "");
                        Summary.WriteWarning(requestingModel, errorMsg);
                        break;
                    default:
                        break;
                }
                return null;
            }
        }

        /// <summary>
        /// Returns the link to the matching resource in the market place if found or creates a new clone copy for future transactions
        /// This allows this action to be performed once to store the link rather than at every transaction
        /// This functionality allows resources not in the market at the start of the simulation to be traded.
        /// </summary>
        /// <param name="resourceType">The resource type to trade</param>
        /// <returns>Whether the search was successful</returns>
        public IResourceWithTransactionType LinkToMarketResourceType(CLEMResourceTypeBase resourceType)
        {
            if (!(this.Parent is Market))
            {
                throw new ApsimXException(this, $"Logic error in code. Trying to link a resource type [r={resourceType.Name}] from the market with the same market./nThis is a coding issue. Please contact the developers");
            }

            // find parent group type
            ResourceBaseWithTransactions parent = (resourceType as Model).Parent as ResourceBaseWithTransactions;
            ResourceBaseWithTransactions resGroup = GetResourceGroupByType(parent.GetType()) as ResourceBaseWithTransactions;
            if (resGroup is null)
            {
                // add warning the market is not currently trading in this resource
                string zoneName = FindAncestor<Zone>().Name;
                string warn = $"[{zoneName}] is currently not accepting resources of type [r={parent.GetType().ToString()}]\r\nOnly resources groups provided in the [r=ResourceHolder] in the simulation tree will be traded.";
                if (!Warnings.Exists(warn) & Summary != null)
                {
                    Summary.WriteWarning(this, warn);
                    Warnings.Add(warn);
                }
                return null;
            }

            // TODO: do some group checks. land units, currency

            // TODO: if market and looking for finance only return or create "Bank"

            // find resource type in group
            object resType = resGroup.FindChild< IResourceWithTransactionType >((resourceType as IModel).Name);
            if (resType is null)
            {
                // clone resource: too many problems with linked events to clone these objects and setup again
                // it will be the responsibility of the user to ensure the resources and details are in the market
                if (resType is null)
                {
                    // add warning the market does not have the resource
                    string zoneName = FindAncestor<Zone>().Name;
                    string warn = $"The resource [r={resourceType.Parent.Name}.{resourceType.Name}] does not exist in [m={this.Parent.Name}].\r\nAdd resource and associated components to the market to permit trading.";
                    if (!Warnings.Exists(warn) & Summary != null)
                    {
                        Summary.WriteWarning(this, warn);
                        Warnings.Add(warn);
                    }
                    return null;
                }
                else
                {
                    (resType as IModel).Parent = resGroup;
                    (resType as CLEMModel).CLEMParentName = resGroup.CLEMParentName;
                    // add new resource type
                    resGroup.AddNewResourceType(resType as IResourceWithTransactionType);
                }
            }
            return resType as IResourceWithTransactionType;
        }

        /// <summary>
        /// Get the Resource Group for Products
        /// </summary>
        /// <returns></returns>
        public ProductStore Products()
        {
            return GetGroupByType(typeof(ProductStore)) as ProductStore;
        }

        /// <summary>
        /// Get the Resource Group for Animal Feed
        /// </summary>
        /// <returns></returns>
        public AnimalFoodStore AnimalFoodStore()
        {
            return GetGroupByType(typeof(AnimalFoodStore)) as AnimalFoodStore;
        }

        /// <summary>
        /// Get the Resource Group for OtherAnimals
        /// </summary>
        /// <returns></returns>
        public OtherAnimals OtherAnimalsStore()
        {
            return GetGroupByType(typeof(OtherAnimals)) as OtherAnimals;
        }

        /// <summary>
        /// Get the Resource Group for FoodStore
        /// </summary>
        /// <returns></returns>
        public HumanFoodStore HumanFoodStore()
        {
            return GetGroupByType(typeof(HumanFoodStore)) as HumanFoodStore;
        }

        /// <summary>
        /// Get the Resource Group for GreenhouseGases
        /// </summary>
        /// <returns></returns>
        public GreenhouseGases GreenhouseGases()
        {
            return GetGroupByType(typeof(GreenhouseGases)) as GreenhouseGases;
        }

        /// <summary>
        /// Get the Resource Group for Labour Family
        /// </summary>
        /// <returns></returns>
        public Labour Labour()
        {
            return GetGroupByType(typeof(Labour)) as Labour;
        }

        /// <summary>
        /// Get the Resource Group for Land
        /// </summary>
        /// <returns></returns>
        public Land Land()
        {
            return GetGroupByType(typeof(Land)) as Land;
        }

        /// <summary>
        /// Get the Resource Group for the GrazeFoodStore
        /// </summary>
        /// <returns></returns>
        public GrazeFoodStore GrazeFoodStore()
        {
            return GetGroupByType(typeof(GrazeFoodStore)) as GrazeFoodStore;
        }

        /// <summary>
        /// Get the Resource Group for Ruminant Herd
        /// </summary>
        /// <returns></returns>
        public RuminantHerd RuminantHerd()
        {
            return GetGroupByType(typeof(RuminantHerd)) as RuminantHerd;
        }

        /// <summary>
        /// Get the Resource Group for Finances
        /// </summary>
        /// <returns></returns>
        public Finance FinanceResource()
        {
            return GetGroupByType(typeof(Finance)) as Finance;
        }

        /// <summary>An event handler to allow us to initialise ourselves.</summary>
        /// <param name="sender">The sender.</param>
        /// <param name="e">The <see cref="EventArgs"/> instance containing the event data.</param>
        [EventSubscribe("Commencing")]
        private void OnSimulationCommencing(object sender, EventArgs e)
        {
            // if this isn't a marketplace try find a shared market
            if(!(this.Parent is Market))
            {
                IModel parentSim = FindAncestor<Simulation>();
                FoundMarket = parentSim.FindAllChildren<Market>().Where(a => a.Enabled).FirstOrDefault();
            }
            else
            {
                FoundMarket = this.Parent as Market;
            }
            InitialiseResourceGroupList();
        }

        /// <summary>
        /// Performs the transmutation of resources into a required resource
        /// </summary>
        public void TransmutateShortfall(List<ResourceRequest> requests, bool queryOnly)
        {
            List<ResourceRequest> shortfallRequests = requests.Where(a => a.Required > a.Available).ToList();

            // Search through all limited resources and determine if transmutation available
            foreach (ResourceRequest request in shortfallRequests)
            {
                // Check if transmutation would be successful 
                if (request.AllowTransmutation && (queryOnly || request.TransmutationPossible))
                {
                    // get resource type
                    IModel model = request.Resource as IModel;
                    if (model is null)
                    {
                        model = this.GetResourceItem(request.ActivityModel, request.ResourceType, request.ResourceTypeName, OnMissingResourceActionTypes.Ignore, OnMissingResourceActionTypes.Ignore) as IModel;
                    }
                    if (model != null)
                    {
                        // get the resource holder to use for this request
                        // not it is either this class or the holder for the market place required.
                        ResourcesHolder resHolder = model.Parent.Parent as ResourcesHolder;

                        // check if transmutations provided
                        foreach (Transmutation trans in model.FindAllChildren<Transmutation>())
                        {
                            double unitsNeeded = 0;
                            // check if resources available for activity and transmutation
                            foreach (ITransmutationCost transcost in trans.FindAllChildren<IModel>().Where(a => a is ITransmutationCost).Cast<ITransmutationCost>())
                            {
                                double unitsize = trans.AmountPerUnitPurchase;
                                if (transcost is TransmutationCostUsePricing)
                                {
                                    // use pricing details if needed
                                    unitsize = (transcost as TransmutationCostUsePricing).Pricing.PacketSize;
                                }
                                unitsNeeded = Math.Ceiling((request.Required - request.Available) / unitsize);

                                double transmutationCost;
                                if (transcost is TransmutationCostUsePricing)
                                {
                                    // use pricing details if needed
                                    transmutationCost = unitsNeeded * (transcost as TransmutationCostUsePricing).Pricing.PricePerPacket;
                                }
                                else
                                {
                                    transmutationCost = unitsNeeded * transcost.CostPerUnit;
                                }

                                // get transcost resource
                                IResourceType transResource = null;
                                if (transcost.ResourceType.Name != "Labour")
                                {
                                    transResource = resHolder.GetResourceItem(request.ActivityModel, transcost.ResourceTypeName, OnMissingResourceActionTypes.Ignore, OnMissingResourceActionTypes.Ignore) as IResourceType;
                                }

                                if (!queryOnly)
                                {
                                    // remove cost
                                    // create new request for this transmutation cost
                                    ResourceRequest transRequest = new ResourceRequest
                                    {
                                        RelatesToResource = request.ResourceTypeName,
                                        Required = transmutationCost,
                                        ResourceType = transcost.ResourceType,
                                        ActivityModel = request.ActivityModel,
                                        Category = "Transmutation",
                                    };

                                    // used to pass request, but this is not the transmutation cost

                                    if (transcost.ResourceType.Name == "Labour")
                                    {
                                        transRequest.ResourceType = typeof(Labour);
                                        transRequest.FilterDetails = (transcost as IModel).FindAllChildren<LabourFilterGroup>().ToList<object>();
                                        CLEMActivityBase.TakeLabour(transRequest, true, transRequest.ActivityModel, this, OnPartialResourcesAvailableActionTypes.UseResourcesAvailable);
                                    }
                                    else
                                    {
                                        transResource.Remove(transRequest);
                                    }
                                }
                                else
                                {
                                    double activityCost = requests.Where(a => a.ResourceType == transcost.ResourceType && a.ResourceTypeName == transcost.ResourceTypeName).Sum(a => a.Required);
                                    if (transmutationCost + activityCost <= transResource.Amount)
                                    {
                                        request.TransmutationPossible = true;
                                        break;
                                    }
                                }
                            }
                            if(!queryOnly)
                            {
                                // Add resource
                                (model as IResourceType).Add(unitsNeeded * trans.AmountPerUnitPurchase, request.ActivityModel, request.ResourceTypeName, "Transmutation");
                            }
                        }
                    }

                }

            }
        }

        #region validation

        /// <summary>
        /// Validate object
        /// </summary>
        /// <param name="validationContext"></param>
        /// <returns></returns>
        public IEnumerable<ValidationResult> Validate(ValidationContext validationContext)
        {
            var results = new List<ValidationResult>();

            var t = this.Children.Where(a => a.GetType().FullName != "Models.Memo").GroupBy(a => a.GetType()).Where(b => b.Count() > 1);

            // check that only one instance of each resource group is present
            foreach (var item in this.Children.Where(a => a.GetType().FullName != "Models.Memo").GroupBy(a => a.GetType()).Where(b => b.Count() > 1))
            {
                string[] memberNames = new string[] { item.Key.FullName };
                results.Add(new ValidationResult(String.Format("Only one (1) instance of any resource group is allowed in the Resources Holder. Multiple Resource Groups [{0}] found!", item.Key.FullName), memberNames));
            }
            return results;
        }


        #endregion

        #region descriptive summary

        /// <summary>
        /// Provides the description of the model settings for summary (GetFullSummary)
        /// </summary>
        /// <param name="formatForParentControl">Use full verbose description</param>
        /// <returns></returns>
        public override string ModelSummary(bool formatForParentControl)
        {
            return "<h1>Resources summary</h1>";
        }

        /// <summary>
        /// Provides the closing html tags for object
        /// </summary>
        /// <returns></returns>
        public override string ModelSummaryOpeningTags(bool formatForParentControl)
        {
            return "\r\n<div class=\"resource\" style=\"opacity: " + SummaryOpacity(formatForParentControl).ToString() + "\">";
        }

        /// <summary>
        /// Provides the closing html tags for object
        /// </summary>
        /// <returns></returns>
        public override string ModelSummaryClosingTags(bool formatForParentControl)
        {
            return "\r\n</div>";
        }

        #endregion
    }
}
