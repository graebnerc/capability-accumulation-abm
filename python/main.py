import sys; sys.path.insert(0, './python')
from model import Model
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
file_marker = "[" + str(os.path.basename(__file__)) + "]: "

class Main:
    """Runs a model for several iterations and computes results
    """

    def __init__(self, parameter_filename, iterations, output_folder):
        print(file_marker + "Parameter file: " + parameter_filename)
        print(file_marker + "Called __init__ of main.py")
        parameters = json.load(open(parameter_filename))
        assert isinstance(parameters, dict), "Parameter file must be dict, \
        not {}".format(type(parameters))
        assert isinstance(iterations, int), "Nb of iterations must be int, \
        not {}".format(type(iterations))

        self.parameters = parameters
        self.base_name = parameter_filename[18:-5]
        self.outcome_filename = output_folder + "/" + self.base_name + ".feather"
        self.outcome_filename_csv = output_folder + "/" + self.base_name + ".csv"
        
        self.results = []
        self.results_dist = []
        self.current_id = 1
        while self.current_id <= iterations:
            print(file_marker + "Start iteration ", self.current_id, " of ", 
                  iterations)
            self.current_model = Model(identifier=self.current_id, 
                                       parameters=self.parameters, 
                                       base_name=self.base_name)
            if self.current_id%10 == 0:
                self.current_model.visualize_product_space(
                    mark_firms=True, 
                    individual_start=output_folder + "/", 
                    individual_end="start")
            self.current_model.run()
            if self.current_id%10 == 0:
                self.current_model.visualize_product_space(
                    mark_firms=True, 
                    individual_start=output_folder + "/", 
                    individual_end="end")
            model_results = self.current_model.return_results()
            self.results.append(model_results["dynamics"])
            self.results_dist.append(model_results["distributions"])
            self.current_id += 1
        self.results_frame = pd.concat(self.results, ignore_index=True)
        self.tmax_frame = self.get_tfinal(self.results_frame)
        self.results_frame = self.aggregate_results(self.results_frame)
        self.results_dist_frame = pd.concat(self.results_dist, 
                                            ignore_index=True)
        self.save_data()
        
    def get_tfinal(self, full_data_frame):
        """Saves final time step for each model run
        """
        final_data_frame = \
            full_data_frame[full_data_frame["t"]==max(full_data_frame["t"])]
        final_data_frame.reset_index(inplace=True)
        return final_data_frame
    
    def aggregate_results(self, full_data_frame):
        """
        Takes the results for individual model runs and aggregates those 
        with the same parameter specification.
        
        Attributes
        ----------
        name : pd.DataFrame
            A data frame containing the results of all model runs.
            Usually created via: `pd.concat(self.results, ignore_index=True)`
            This should be run automatically one line before calling this 
            function.

        Returns
        -------
        pd.DataFrame
            Data frame with means and standard deviations for all 
            relevant state variables of the model.
        """
        grouping_parameters = ["t", "n_firms", "n_workers", "n_products",
                               "prod_space_struc", "delta_coefficient",
                               "prod_sp_val_alct", "prod_sp_val_dist",
                               "initial_capital_stock", "initial_capabilities",
                               "nominal_demand", "fin_regime", "cap_prod",
                               "max_rd_success", "depreciation_rate", "max_info"]
        full_data_frame_agg = full_data_frame.groupby(grouping_parameters).agg(
            output_all_prdcts_mean=("output_all_prdcts", np.mean),
            output_all_prdcts_std=("output_all_prdcts", np.std),
            share_produced_products_mean=("share_produced_products", np.mean),
            share_produced_products_std=("share_produced_products", np.std),
            comp_produced_prod_mean_mean=("comp_produced_prod_mean", np.mean),
            comp_produced_prod_mean_std=("comp_produced_prod_mean", np.std),
            price_produced_prod_mean_mean=("price_produced_prod_mean", np.mean),
            price_produced_prod_mean_std=("price_produced_prod_mean", np.std),
            price_produced_prod_std_mean=("price_produced_prod_std", np.mean),
            price_produced_prod_std_std=("price_produced_prod_std", np.std),
            correlation_price_complexity_mean=(
                "correlation_price_complexity", np.mean),
            correlation_price_complexity_std=(
                "correlation_price_complexity", np.std),
            agg_capital_stock_mean=("agg_capital_stock", np.mean),
            agg_capital_stock_std=("agg_capital_stock", np.std),
            firm_inno_cap_mean_mean=("firm_inno_cap_mean", np.mean),
            firm_inno_cap_mean_std=("firm_inno_cap_mean", np.std),
            firm_inno_cap_sd_mean=("firm_inno_cap_sd", np.mean),
            firm_inno_cap_sd_std=("firm_inno_cap_sd", np.std),
            firm_inno_cap_total_mean=("firm_inno_cap_total", np.mean),
            firm_inno_cap_total_std=("firm_inno_cap_total", np.std),
            firm_spillover_cap_mean_mean=("firm_spillover_cap_mean", np.mean),
            firm_spillover_cap_mean_std=("firm_spillover_cap_mean", np.std),
            firm_spillover_cap_sd_mean=("firm_spillover_cap_sd", np.mean),
            firm_spillover_cap_sd_std=("firm_spillover_cap_sd", np.std),
            firm_spillover_cap_total_mean=("firm_spillover_cap_total", np.mean),
            firm_spillover_cap_total_std=("firm_spillover_cap_total", np.std),
            firm_rd_invest_mean_mean=("firm_rd_invest_mean", np.mean),
            firm_rd_invest_mean_std=("firm_rd_invest_mean", np.std),
            firm_rd_invest_sd_mean=("firm_rd_invest_sd", np.mean),
            firm_rd_invest_sd_std=("firm_rd_invest_sd", np.std),
            firm_rd_invest_total_mean=("firm_rd_invest_total", np.mean),
            firm_rd_invest_total_std=("firm_rd_invest_total", np.std),
            firm_ac_invest_mean_mean=("firm_ac_invest_mean", np.mean),
            firm_ac_invest_mean_std=("firm_ac_invest_mean", np.std),
            firm_ac_invest_sd_mean=("firm_ac_invest_sd", np.mean),
            firm_ac_invest_sd_std=("firm_ac_invest_sd", np.std),
            firm_ac_invest_total_mean=("firm_ac_invest_total", np.mean),
            firm_ac_invest_total_std=("firm_ac_invest_total", np.std),
            firm_k_invest_mean_mean=("firm_k_invest_mean", np.mean),
            firm_k_invest_mean_std=("firm_k_invest_mean", np.std),
            firm_k_invest_sd_mean=("firm_k_invest_sd", np.mean),
            firm_k_invest_sd_std=("firm_k_invest_sd", np.std),
            firm_k_invest_total_mean=("firm_k_invest_total", np.mean),
            firm_k_invest_total_std=("firm_k_invest_total", np.std),
            firm_profits_mean_mean=("firm_profits_mean", np.mean),
            firm_profits_mean_std=("firm_profits_mean", np.std),
            firm_profits_sd_mean=("firm_profits_sd", np.mean),
            firm_profits_sd_std=("firm_profits_sd", np.std),
            firm_profits_total_mean=("firm_profits_total", np.mean),
            firm_profits_total_std=("firm_profits_total", np.std),
            firm_gross_revenues_mean_mean=("firm_gross_revenues_mean", np.mean),
            firm_gross_revenues_mean_std=("firm_gross_revenues_mean", np.std),
            firm_gross_revenues_sd_mean=("firm_gross_revenues_sd", np.mean),
            firm_gross_revenues_sd_std=("firm_gross_revenues_sd", np.std),
            firm_gross_revenues_total_mean=("firm_gross_revenues_total", np.mean),
            firm_gross_revenues_total_std=("firm_gross_revenues_total", np.std),
            firm_accounts_mean_mean=("firm_accounts_mean", np.mean),
            firm_accounts_mean_std=("firm_accounts_mean", np.std),
            firm_accounts_sd_mean=("firm_accounts_sd", np.mean),
            firm_accounts_sd_std=("firm_accounts_sd", np.std),
            firm_accounts_total_mean=("firm_accounts_total", np.mean),
            firm_accounts_total_std=("firm_accounts_total", np.std),
            firm_share_waiting_mean=("firm_share_waiting", np.mean),
            firm_share_waiting_std=("firm_share_waiting", np.std),
            firm_share_remained_mean=("firm_share_remained", np.mean),
            firm_share_remained_std=("firm_share_remained", np.std),
            firm_deposits_mean_mean=("firm_deposits_mean", np.mean),
            firm_deposits_mean_std=("firm_deposits_mean", np.std),
            firm_assets_mean_mean=("firm_assets_mean", np.mean),
            firm_assets_sd_mean=("firm_assets_sd", np.mean),
            firm_assets_total_mean=("firm_assets_total", np.mean),
            firm_liabilities_mean_mean=("firm_liabilities_mean", np.mean),
            firm_liabilities_sd_mean=("firm_liabilities_sd", np.mean),
            firm_liabilities_total_mean=("firm_liabilities_total", np.mean),
            firm_net_worth_mean_mean=("firm_net_worth_mean", np.mean),
            firm_net_worth_sd_mean=("firm_net_worth_sd", np.mean),
            firm_net_worth_total_mean=("firm_net_worth_total", np.mean),
            bank_assets_mean_mean=("bank_assets_mean", np.mean),
            bank_assets_sd_mean=("bank_assets_sd", np.mean),
            bank_assets_total_mean=("bank_assets_total", np.mean),
            bank_liabilities_mean_mean=("bank_liabilities_mean", np.mean),
            bank_liabilities_sd_mean=("bank_liabilities_sd", np.mean),
            bank_liabilities_total_mean=("bank_liabilities_total", np.mean),
            bank_net_worth_mean_mean=("bank_net_worth_mean", np.mean),
            bank_net_worth_sd_mean=("bank_net_worth_sd", np.mean),
            bank_net_worth_total_mean=("bank_net_worth_total", np.mean),
            firms_final_spot_mean=("firms_final_spot", np.mean),
            firms_share_deposits_invested_mean_mean=(
                "firm_share_investment_mean", np.mean),
            firms_share_deposits_invested_sd_mean=(
                "firm_share_investment_sd", np.mean),
            firm_assets_mean_std=("firm_assets_mean", np.std),
            firm_assets_sd_std=("firm_assets_sd", np.std),
            firm_assets_total_std=("firm_assets_total", np.std),
            firm_liabilities_mean_std=("firm_liabilities_mean", np.std),
            firm_liabilities_sd_std=("firm_liabilities_sd", np.std),
            firm_liabilities_total_std=("firm_liabilities_total", np.std),
            firm_net_worth_mean_std=("firm_net_worth_mean", np.std),
            firm_net_worth_sd_std=("firm_net_worth_sd", np.std),
            firm_net_worth_total_std=("firm_net_worth_total", np.std),
            bank_assets_mean_std=("bank_assets_mean", np.std),
            bank_assets_sd_std=("bank_assets_sd", np.std),
            bank_assets_total_std=("bank_assets_total", np.std),
            bank_liabilities_mean_std=("bank_liabilities_mean", np.std),
            bank_liabilities_sd_std=("bank_liabilities_sd", np.std),
            bank_liabilities_total_std=("bank_liabilities_total", np.std),
            bank_net_worth_mean_std=("bank_net_worth_mean", np.std),
            bank_net_worth_sd_std=("bank_net_worth_sd", np.std),
            bank_net_worth_total_std=("bank_net_worth_total", np.std),
            firms_final_spot_std=("firms_final_spot", np.std),
            firms_share_deposits_invested_mean_std=(
                "firm_share_investment_mean", np.std),
            firms_share_deposits_invested_sd_std=(
                "firm_share_investment_sd", np.std)
        )
        full_data_frame_agg.reset_index(inplace=True)
        return full_data_frame_agg        

    def save_data(self):
        """Saves results in a csv file
        """
        print(file_marker + "Start saving data...", end="")
        dist_outcome_filename = \
            self.outcome_filename.replace(".feather", "_dist.feather")
        tmax_outcome_filename = \
            self.outcome_filename.replace(".feather", "_tmax.feather")
        tmax_outcome_filename_csv = \
            self.outcome_filename.replace(".feather", "_tmax.csv")
        self.results_frame.reset_index(drop=True).to_feather(self.outcome_filename)
        self.results_frame.reset_index(drop=True).to_csv(self.outcome_filename_csv)
        # self.results_dist_frame.reset_index(drop=True).to_feather(dist_outcome_filename)
        self.tmax_frame.reset_index(drop=True).to_feather(tmax_outcome_filename)
        self.tmax_frame.to_csv(tmax_outcome_filename_csv)
        print(file_marker + "complete!")
        print(file_marker + "Outcome saved in: {}".format(self.outcome_filename))
        print(file_marker + "Final-t outcome saved in: {}".format(
            tmax_outcome_filename))

    @staticmethod
    def get_colors(cm_name, nb_cols):
        """Gives a list of color codes from a given color map.

        Parameters
        ----------
        cm_name : str
            Name of a color map used in matplotlib.pyplot
        nb_cols : int
            Number of different colors to be returned.

        Returns
        -------
        list
            A list with `nb_cols` color codes from color map `cm_name`.
        """
        cmap_object = plt.get_cmap(cm_name)
        col_codes = cmap_object(np.linspace(1, 256, nb_cols) / 100)
        return col_codes
    