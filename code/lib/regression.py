from .correction_functions import clean_IRENA_summary
from .util import *


def pyomo_regression_model():
    """
    This function returns an abstract pyomo model of a constrained least square problem for time series fitting to
    match model FLHs and minimize difference error with model time series.
    
    :return model: Abstract pyomo model.
    :rtype: pyomo object
    """
    model = pyo.AbstractModel()
    model.s = pyo.Set()
    model.q = pyo.Set()
    model.t = pyo.Set()

    model.FLH = pyo.Param()
    model.shape = pyo.Param(model.t)

    model.TS = pyo.Param(model.s, model.q, model.t)
    model.coef = pyo.Var(model.s, model.q, domain=pyo.NonNegativeReals)

    def constraint_FLH(model):
        FLH = 0
        for s in model.s:
            for q in model.q:
                tempTS = 0
                for t in model.t:
                    tempTS = tempTS + model.TS[s, q, t]
                FLH = FLH + pyo.prod([model.coef[s, q], tempTS])
        return FLH == model.FLH

    def constraint_sum(model):
        sum = 0
        for s in model.s:
            for q in model.q:
                sum = sum + model.coef[s, q]
        return sum == 1

    def obj_expression(model):
        Error = 0
        for s in model.s:
            for q in model.q:
                for t in model.t:
                    Error = Error + (pyo.prod([model.coef[s, q], model.TS[s, q, t]]) - model.shape[t]) ** 2
        return Error

    model.OBJ = pyo.Objective(rule=obj_expression)
    model.constraint_FLH = pyo.Constraint(rule=constraint_FLH)
    model.constraint_sum = pyo.Constraint(rule=constraint_sum)
    return model


def clean_FLH_regression(paths, param):
    """
    This function creates a CSV file containing the model FLH used for regression. If the region is present in the
    IRENA database, then the FLH are extracted directly from there. In case it is not present, a place holder for the
    regions is written in the csv file and it is the user's responsibility to fill in an appropriate value.
    The function will warn the user, and print all regions that are left blank.

    :param param: Dictionary of dictionaries containing the list of regions.
    :type param: dict
    :param paths: Dictionary of dictionaries containing the paths to *IRENA_summary*, *IRENA_dict*.
    :type paths: dict

    :return missing: List of string of the missing regions. The CSV file for the the FLH needed for the regression is saved directly in
        the given path, along with the corresponding metadata in a JSON file.
    :rtype: list of str
    :raise Missing Regions: No FLH values exist for certain regions.
    """
    # Read IRENA summary
    if not os.path.isfile(paths["IRENA_summary"]):
        clean_IRENA_summary(paths, param)
    IRENA_summary = pd.read_csv(paths["IRENA_summary"], sep=";", decimal=",", index_col=[0, 1])

    # Load IRENA dictionary
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";")
    IRENA_dict.dropna(inplace=True)
    IRENA_dict.set_index(["NAME_SHORT"], inplace=True)
    IRENA_dict = IRENA_dict["Countries shapefile"].to_dict()

    # Setup FLH_regression dataframe
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    FLH_regression = pd.DataFrame(columns=["WindOn", "WindOff", "PV", "CSP"], index=list_regions, dtype=float)

    # Fill in FLH_regression dataframe
    missing = []
    for country in list_regions:
        # Country/region found in IRENA_dict and IRENA_summary
        if country in IRENA_dict.keys():
            FLH_regression.loc[country] = [
                IRENA_summary.loc[(IRENA_dict[country], "Onshore wind energy"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Offshore wind energy"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Solar photovoltaic"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Concentrated solar power"), "FLH (h)"],
            ]
        # Missing country/region, require user input
        else:
            missing.append(country)

    # Save FLH_regression
    FLH_regression.to_csv(paths["FLH_regression"], sep=";", decimal=",", index=True)
    create_json(
        paths["FLH_regression"],
        param,
        ["author", "comment", "region_name", "subregions_name", "year", "Crd_all"],
        paths,
        ["IRENA_dict", "IRENA_summary"],
    )
    print("files saved: " + paths["FLH_regression"])

    # Return missing countries/regions
    warn(
        "The following countries/regions are not present in the IRENA Database: "
        + ",".join(missing)
        + ".\nTheir corresponding FLH have been left blank.",
        UserWarning,
    )
    return missing


def clean_TS_regression(paths, param, tech):
    """
    This function creates a CSV file containing the model time series used for regression. If the region is present in
    the EMHIRES text files then the TS is extracted directly from it. If the region is not present in the EMHIRES text
    files, the highest FLH generated TS is used instead and is scaled to match IRENA FLH if the IRENA FLH are available.

    :param paths: Dictionary containing paths to EMHIRES text files.
    :type paths: dict
    :param param: Dictionary containing the *FLH_regression* dataframe, list of subregions contained in shapefile, and year.
    :type param: dict

    :return: The time series used for the regression are saved directly in the given path, along with the corresponding metadata in a JSON file. 
    :rtype: None
    :raise Missing FLH: FLH values are missing for at least one region. No scaling is applied to the time series for those regions.
    :raise Missing EMHIRES: EMHIRES database is missing, generated timeseries will be used as model for all regions.
    """

    # load IRENA FLH data
    irena = pd.read_csv(paths["FLH_regression"], sep=";", decimal=",", index_col=0)

    # Find intersection between desired regions and irena regions
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    list_regions = sorted(list(set(list_regions).intersection(set(irena.index))))

    # Create TS_regression dataframe
    TS_regression = pd.DataFrame(index=range(1, 8761), columns=list_regions)
    if os.path.isfile(paths[tech]["EMHIRES"]):
        # Load EMHIRES data for desired year
        if tech in ["PV", "CSP"]:
            date_index = pd.date_range(start="1/1/1986", end="1/1/2016", freq="H", closed="left")
            EMHIRES = pd.read_csv(paths[tech]["EMHIRES"], " ")
            EMHIRES = EMHIRES.set_index(date_index)
            EMHIRES = EMHIRES.loc["1/1/" + str(param["year"]) : "1/1/" + str(param["year"] + 1)]
        else:
            EMHIRES = pd.read_csv(paths[tech]["EMHIRES"], "\t")
            EMHIRES = EMHIRES[EMHIRES["Year"] == param["year"]].reset_index()
            EMHIRES = EMHIRES.drop(["index", "Time step", "Date", "Year", "Month", "Day", "Hour"], axis=1)

        # Find intersection between EMHIRES and list_regions
        intersect_regions = sorted(list((set(list_regions).intersection(set(EMHIRES.columns)))))
    else:
        intersect_regions = []
        warn("WARNING: EMHIRES database unavailable for technology " + tech, UserWarning)
    # Load setting combinations
    settings = combinations_for_regression(paths, param, tech)
    if not settings[0]:
        settings = settings[1]
    else:
        settings = settings[0]
    nanval = ""
    for region in list_regions:
        IRENA_FLH = irena.loc[region, tech]
        # Region is present in both EMHIRES and IRENA
        if region in intersect_regions:
            # Scale EMHIRES TS to IRENA
            if not np.isnan(IRENA_FLH):
                TS_regression[region] = (EMHIRES[region] * (IRENA_FLH / sum(EMHIRES[region]))).values
            else:
                TS_regression = EMHIRES[region].values
                nanval = nanval + region + ", "
        # Region is not present in EMHIRES, use scaled generated TS instead
        else:
            # Load generated TS and scale it with IRENA FLH
            GenTS = read_generated_TS(paths, param, tech, settings, region)
            # Find highest FLH valued TS
            settings_sorted = np.array(
                pd.DataFrame((np.nansum(GenTS[key]) for key in GenTS.keys()), index=settings, columns=["FLH_all_quant"])
                .sort_values(by="FLH_all_quant", ascending=0)
                .index
            )
            GenTS["TS_Max"] = GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))]
            # Scale max TS to IRENA
            if not np.isnan(IRENA_FLH):
                TS_regression[region] = (GenTS["TS_Max"] * (IRENA_FLH / GenTS["TS_Max"].sum())).values
            else:
                TS_regression[region] = GenTS["TS_Max"].values
                nanval = nanval + region + ", "

    # Save TS_regression as CSV
    TS_regression.to_csv(paths[tech]["TS_regression"], sep=";", decimal=",", index=True)
    create_json(
        paths[tech]["TS_regression"], param, ["author", "comment", tech, "region_name", "subregions_name", "year"], paths, ["FLH_regression", tech]
    )
    if nanval != "":
        warn("Missing FLH for regions: " + nanval.rstrip(", ") + "\nTime series have not been scaled for these regions")
    print("files saved: " + paths[tech]["TS_regression"])


def get_regression_coefficients(paths, param, tech):
    """
    This function solves the following optimization problem: A combination of quantiles, hub heights or orientations is to be found, so that
    the error to a given historical time series (e.g. from EMHIRES for European countries) is minimized, while
    constraining the FLH to match a given value (for example from IRENA). The settings of the combinations can be
    defined by the user.

    The function starts by identifying the existing settings (hub heights, orientations) and quantiles.
    If the combinations of time series requested by the user cannot be found, a warning is raised.

    It later runs the optimization and identifies the subregions for which a solution was found. If the optimization
    is infeasible (too high or too low FLH values compared to the reference to be matched), the time series with the closest
    FLH to the reference value is used in the final output.

    The output consists of coefficients between 0 and 1 that could be multiplied later with the individual time series
    in :mod:`time_series.generate_stratified_timeseries`. The sum of the coefficients for each combination is equal to 1.

    :param paths: Dictionary including the paths to the time series for each subregion, technology setting, and quantile, to the output paths for the coefficients.
    :type paths: dict
    :param param: Dictionary including the dictionary of regression parameters, quantiles, and year.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return:
        The regression parameters (e.g. IRENA FLH and EMHIRES TS) are copied under *regression_in* folder, and
        the regression coefficients are saved in a CSV file under *regression_out* folder, along with the metadata
        in a JSON file.
    :rtype: None
    :raise Missing Data: No time series present for technology *tech*.
    :raise Missing Data for Setting: Missing time series for desired settings (hub heights / orientations).
    """
    timecheck("Start")
    year = str(param["year"])

    try:
        combinations = combinations_for_regression(paths, param, tech)
        param["combinations"] = combinations
    except UserWarning:
        timecheck("End")
        return
    if combinations is None:
        timecheck("End")
        return

    # Display the combinations of settings to be used
    if tech in ["WindOn", "WindOff"]:
        print("Combinations of hub heights to be used for the regression: ", combinations)
    elif tech in ["PV"]:
        print("Orientations to be used for the regression: ", combinations)

    # Create FLH file for regression
    if not os.path.isfile(paths["FLH_regression"]):
        clean_FLH_regression(paths, param)

    # Create TS file for regression
    if not os.path.isfile(paths[tech]["TS_regression"]):
        clean_TS_regression(paths, param, tech)

    FLH, TS_reg = check_regression_model(paths, tech)

    param["FLH_regression"] = FLH
    param["TS_regression"] = TS_reg

    # Find intersection between FLH and shapefile subregions
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    list_regions = sorted(list(set(list_regions).intersection(set(FLH.index))))

    # loop over all combinations and regions
    for settings in combinations:
        # Summary variables
        summary = None
        nodata = ""
        no_sol_high = ""
        no_sol_low = ""
        solution = ""
        status = 0
        print("Regions under study : ", list_regions)
        for reg in list_regions:
            # Show progress of the simulation
            status = status + 1
            display_progress("Regression Coefficients " + tech + " " + param["subregions_name"], (len(list_regions), status))

            region_data = regmodel_load_data(paths, param, tech, settings, reg)

            # Skip regions not present in the generated TS
            if region_data is None:
                nodata = nodata + reg + ", "
                settings_sorted = settings
                continue

            settings_sorted = region_data[None]["s"][None].tolist()

            if region_data[None]["IRENA_best_worst"] == (True, True):

                # create model instance
                solver = SolverFactory(param["regression"]["solver"])
                model = pyomo_regression_model()
                regression = model.create_instance(region_data)

                # solve model and return results
                solver.solve(regression)

                # Retrieve results
                r = np.zeros((len(param["quantiles"]), len(settings_sorted)))
                c = 0
                for q in param["quantiles"]:
                    p = 0
                    for s in settings_sorted:
                        r[c, p] = pyo.value(regression.coef[s, q])
                        p += 1
                    c += 1
                r[r < 10 ** (-5)] = 0
                solution = solution + reg + ", "

            elif region_data[None]["IRENA_best_worst"] == (False, True):
                # Select best TS (highest height, highest quantile)
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)
                r[0, 0] = 1
                no_sol_high = no_sol_high + reg + ", "

            elif region_data[None]["IRENA_best_worst"] == (True, False):
                # Select worst TS (lowest height, lowest quantile)
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)
                r[-1, -1] = 1
                no_sol_low = no_sol_low + reg + ", "

            else:
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)

            if settings_sorted != [0]:
                result = pd.DataFrame(r, param["quantiles"], (reg + "_" + str(s) for s in settings_sorted))
            else:
                result = pd.DataFrame(r, param["quantiles"], [reg])

            if summary is None:
                summary = result
            else:
                summary = pd.concat([summary, result], axis=1)

        # Print Regression Summary
        if solution != "":
            print("\nA solution was found for the following regions: " + solution.rstrip(", "))
        if no_sol_low != "":
            print("\nNo Solution was found for the following regions because they are too high: " + no_sol_low.rstrip(", "))
        if no_sol_high != "":
            print("\nNo Solution was found for the following regions because they are too low: " + no_sol_high.rstrip(", "))
        if nodata != "":
            print("\nNo data was available for the following regions: " + nodata.rstrip(", "))

        if summary is None:
            return
        st = ""
        for setting in settings_sorted:
            st = st + str(setting) + "_"

        summary.to_csv(paths[tech]["Regression_coefficients"] + st + year + ".csv", sep=";", decimal=",")
        create_json(
            paths[tech]["Regression_coefficients"],
            param,
            ["author", "comment", tech, "region_name", "subregions_name", "quantiles", "regression", "year", "Crd_all"],
            paths,
            ["spatial_scope", "subregions"],
        )
        print("\nfiles saved: " + paths[tech]["Regression_coefficients"] + st + year + ".csv")

    timecheck("End")


def read_generated_TS(paths, param, tech, settings, subregion):
    """
    This function returns a dictionary containing the available time series generated by the script based on
    the desired technology and settings.

    :param paths: Dictionary including output folder for regional analysis.
    :type paths: dict
    :param param: Dictionary including list of subregions and year.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str
    :param settings: List of lists containing setting combinations.
    :type settigns: list
    :param subregion: Name of the subregion.
    :type subregion: str

    :return GenTS: Dictionary of time series indexed by setting and quantile.
    :rtype: dict
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    bef_setting = paths["regional_analysis"] + subregions + "_" + tech + "_"
    aft_setting = "_TS_" + year + ".csv"

    # Setup the data dictionary for generated TS for each quantile
    GenTS = {}

    for setting in settings:
        TS_Temp = pd.read_csv(bef_setting + str(setting) + aft_setting, sep=";", decimal=",", dtype=str)

        filter_reg = [col for col in TS_Temp if col.startswith(subregion)]
        # Remove undesired regions
        TS_Temp = TS_Temp[filter_reg]

        # Exit function if subregion is not present in TS files
        if TS_Temp.empty:
            return None

        TS_Temp.columns = TS_Temp.iloc[0]
        TS_Temp = TS_Temp.drop(0)
        # Replace ',' with '.' for float conversion
        for q in range(0, len(TS_Temp.columns)):
            TS_Temp.iloc[:, q] = TS_Temp.iloc[:, q].apply(lambda x: x.replace(",", ".")).astype(float)
            TS_Temp.columns.name = ""
            TS_Temp.reset_index(inplace=True, drop=True)
        GenTS[str(setting)] = TS_Temp

    return GenTS


def regmodel_load_data(paths, param, tech, settings, subregion):
    """
    This function returns a dictionary used to initialize a pyomo abstract model for the regression analysis
    of each region.

    :param paths: Dictionary of dictionaries containing the paths to the CSV time series files.
    :type paths: dict
    :param param: Dictionary of dictionaries contating IRENA's region list, FLHs and EMHIRES model timeseries.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str
    :param settings: List of all the settings (hub heights/orientations) to be used in the regression.
    :type settings: list
    :param subregion: Name of subregion.
    :type subregion: str

    :return data: Dictionary containing regression parameters.
    :rtype: dict
    """
    subregions = param["subregions_name"]
    year = str(param["year"])
    FLH = param["FLH_regression"]
    time = range(1, 8761)

    # Setup dataframe for IRENA
    model_FLH = FLH.loc[subregion, tech]
    if np.isnan(model_FLH):
        return None

    # Read data from output folder
    GenTS = read_generated_TS(paths, param, tech, settings, subregion)
    if GenTS is None:
        return None

    # reorder hubheights to go from max TS to min TS:
    settings_sorted = np.array(
        pd.DataFrame((np.nansum(GenTS[key]) for key in GenTS.keys()), index=settings, columns=["FLH_all_quant"])
        .sort_values(by="FLH_all_quant", ascending=0)
        .index
    )
    GenTS["TS_Max"] = np.nansum(GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))])
    GenTS["TS_Min"] = np.nansum(GenTS[str(settings_sorted[-1])]["q" + str(np.min(param["quantiles"]))])

    # Check if solution exists
    solution_check = (GenTS["TS_Max"] > model_FLH, GenTS["TS_Min"] < model_FLH)

    # Prepare Timeseries dictionary indexing by height and quantile

    if solution_check == (False, True):
        Timeseries = GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))]

    elif solution_check == (True, False):
        Timeseries = GenTS[str(settings_sorted[-1])]["q" + str(np.min(param["quantiles"]))]

    elif solution_check == (True, True):
        Timeseries = {}

        for s in settings_sorted:
            for q in param["quantiles"]:
                for t in time:
                    Timeseries[(s, q, t)] = np.array(GenTS[str(s)]["q" + str(q)])[t - 1]

    # Setup dataframe for TS Models
    TS_reg = param["TS_regression"]
    ts = np.array(TS_reg[subregion].values)

    TS = {}
    for t in time:
        TS[(t,)] = ts[t - 1]
    # Create data_input dictionary
    data = {
        None: {
            "s": {None: settings_sorted},
            "q": {None: param["quantiles"]},
            "FLH": {None: model_FLH},
            "shape": TS,
            "t": {None: np.array(time)},
            "TS": Timeseries,
            "IRENA_best_worst": solution_check,
            "GenTS": GenTS,
        }
    }
    return data


def combinations_for_regression(paths, param, tech):
    """
    This function reads the list of generated time series for different hub heights and orientations, compares it to the
    user-defined combinations and returns a list of lists containing all the available combinations. The function will
    return a warning if the user input and the available time series are not congruent.

    :param paths: Dictionary of dictionaries containing the paths to the regional analysis output folder.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the subregions name, year, and user-defined combinations.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return combinations: List of combinations for regression.
    :rtype: list
    :raise missing data: If no time series are available for this technology, a warning is raised.
    :raise missing combination: If a hub height or orientation is missing based on user-defined combinations, a warning is raised.
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    # Reads the files present in input folder
    inputfiles = glob(paths["regional_analysis"] + subregions + "_" + tech + "*_TS_" + year + ".csv")

    # Case 1: no files existing
    if len(inputfiles) == 0:
        warn("Generate time series first, before doing the regression!", UserWarning)
        return

    # Get existing settings
    settings_existing = []
    for filename in inputfiles:
        bef_setting = paths["regional_analysis"] + subregions + "_" + tech + "_"
        aft_setting = "_TS_" + year + ".csv"
        settings_existing = settings_existing + [int(filename.replace(bef_setting, "").replace(aft_setting, ""))]
    settings_existing = set(settings_existing)
    print("\nFor technology " + tech + ", time series for the following settings have been detected: ", settings_existing)

    # Get required settings
    combinations = param["regression"][tech].values()
    combinations_sorted = []
    for combi in combinations:
        combinations_sorted = combinations_sorted + [sorted(combi)]
    combinations = combinations_sorted
    settings_required = set([item for sublist in combinations for item in sublist])

    # Case 2: some files are missing
    if not settings_required.issubset(settings_existing):
        print("\nFor technology " + tech + ", time series for the following settings are required: ", settings_required)
        warn("Not all time series are available! Generate the missing time series first, then do the regression.", UserWarning)
        return
    # Case 3: all the required files exist (and nothing more)
    if settings_existing.issubset(settings_required):
        if [] in combinations:
            full_combi = sorted(list(settings_existing))
            if not (full_combi in combinations):
                combinations = combinations + [full_combi]
            combinations.remove([])
    # Case 4: more files exist than those required
    else:
        if [] in combinations:
            full_combi = sorted(list(settings_existing))
            if not (full_combi in combinations):
                combinations = combinations + [full_combi]
            combinations.remove([])
    return combinations


def check_regression_model(paths, tech):
    """
    This function checks the regression model parameters for nan values, and returns the FLH and TS model dataframes.
    If missing values are present in the input CSV files, the users are prompted if they wish to continue or can modify
    the corresponding files.

    :param paths: Dictionary of dictionaries containing the paths to the FLH and TS model regression CSV files.
    :type paths: dict
    :param tech: Technology under study.
    :type tech: str

    :return (FLH, TS_reg): Tuple of pandas dataframes for FLH and TS.
    :rtype: Tuple of pandas dataframes
    """
    while True:
        # Load IRENA data and regions
        FLH = pd.read_csv(paths["FLH_regression"], sep=";", decimal=",", index_col=0)

        # load TS regression file
        TS_reg = pd.read_csv(paths[tech]["TS_regression"], sep=";", decimal=",", index_col=0, header=0)

        # Create filter for nan and 0 values for FLH_regression
        filter_FLH = np.logical_or(np.isnan(FLH[tech]), FLH[tech] == 0)
        reg_nan_null = list(FLH.loc[filter_FLH].index)

        if len(reg_nan_null) != 0:
            print("Missing data:" + ",".join(reg_nan_null))
            ans = input("Some regions are missing FLH data for the technology of choice. Continue ? [y]/n ")
            if ans in ["", "y", "[y]"]:
                # Load IRENA data and regions
                FLH = pd.read_csv(paths["FLH_regression"], sep=";", decimal=",", index_col=0)
                break
        else:
            break

    return FLH, TS_reg
