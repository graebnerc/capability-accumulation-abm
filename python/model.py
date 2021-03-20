import os
import sys
import numpy as np
import pandas as pd
import pdb
import networkx as nx
import matplotlib.pyplot as plt
import math
from collections import Counter
from scipy.stats.stats import pearsonr
from operator import itemgetter
from firm import Firm
from bank import Bank

sys.path.insert(0, './python')
file_marker = "[" + str(os.path.basename(__file__)) + "]: "


class Model:
    def __init__(self, identifier, parameters, base_name):
        """Initializes the model.

        1. Sets relevant parameters as model properties.
            This is not necessary but makes the reporting of intermediate
            results easier during debugging.
        2. Sets up history lists with initial values.
        3. Creates the product space and allocates complexity values.
        4. Creates firms and positions them on the product space.
        5. Creates banks and tests the consistency of model assumptions.
        
        Parameters
        ----------
        identifier : int
            An identifier for the model instance. Will be created within 
            `main.Main`
        parameters : dict
            A dictionary with parameters. Will be created within 
            `main.Main` from the relevant json file
        base_name : str
            The name of the model used for saving results. Will be 
            created in `main.Main` automatically from the name of the 
            json file containing the parameters
        """
        print(file_marker + "Initializing model " + str(identifier) + "...")
        self.id = identifier
        self.base_name = base_name
        self.t = -1  # Time starts in period -1, when product space is created.

        # Set general parameters
        self.parameters = parameters
        self.nb_of_timesteps = parameters["number_of_timesteps"]
        self.n_firms = parameters["number_of_firms"]
        self.n_banks = parameters["number_of_banks"]
        self.n_workers = parameters["number_of_workers"]
        self.delta_coefficient = parameters["delta_coefficient"]
        self.initial_capabilities = parameters['initial_capabilities_equal']
        self.initial_capital_stock = parameters['initial_capital_stock']
        self.nominal_demand = parameters['nominal_demand']
        self.financial_regime_parameter = parameters[
            'financial_regime_parameter']
        self.capital_productivity = parameters["capital_productivity"]
        self.max_rd_success = parameters["max_rd_success"]
        self.depreciation_rate = parameters["depreciation_rate"]
        self.p_innovation = parameters["p_innovation"]
        self.p_spillovers = parameters["p_spillovers"]
        self.p_information = parameters["p_information"]
        self.rate_on_deposits = parameters["rate_on_deposits"]
        self.rate_on_loans = parameters["rate_on_loans"]
        self.payback_rate = parameters["payback_rate"]
        self.share_firms_waiting = [np.nan]
        self.share_firms_remained = [np.nan]

        # Set initial values
        # This is the initial production value for all products
        self.aggregate_output_initial = \
            self.n_firms * self.initial_capabilities \
            * self.initial_capital_stock

        # Set history lists for state variables
        self.agg_capital_stock = [self.initial_capital_stock * self.n_firms]
        self.agg_output_all_prdcts = [self.aggregate_output_initial]
        self.products_produced_share = [np.nan]
        self.products_produced_avg_comp = [np.nan]
        self.firm_size_mean = [np.nan]
        self.firm_size_sd = [np.nan]
        self.prize_avg = [np.nan]
        self.prize_sd = [np.nan]
        self.prize_compl_cor = [np.nan]

        self.firm_innovation_capabilities_mean = [np.nan]
        self.firm_innovation_capabilities_sd = [np.nan]
        self.firm_innovation_capabilities_total = [np.nan]

        self.firm_spillover_capabilities_mean = [np.nan]
        self.firm_spillover_capabilities_sd = [np.nan]
        self.firm_spillover_capabilities_total = [np.nan]

        self.firm_rd_investment_mean = [np.nan]
        self.firm_rd_investment_sd = [np.nan]
        self.firm_rd_investment_total = [np.nan]

        self.firm_k_investment_mean = [np.nan]
        self.firm_k_investment_sd = [np.nan]
        self.firm_k_investment_total = [np.nan]

        self.firm_ac_investment_mean = [np.nan]
        self.firm_ac_investment_sd = [np.nan]
        self.firm_ac_investment_total = [np.nan]

        self.firm_profits_mean = [np.nan]
        self.firm_profits_sd = [np.nan]
        self.firm_profits_total = [np.nan]

        self.firm_gross_revenues_mean = [np.nan]
        self.firm_gross_revenues_sd = [np.nan]
        self.firm_gross_revenues_total = [np.nan]

        self.firm_accounts_mean = [np.nan]
        self.firm_accounts_sd = [np.nan]
        self.firm_accounts_total = [np.nan]

        self.firm_profit_rate_mean = [np.nan]
        self.firm_profit_rate_sd = [np.nan]
        self.firm_profit_rate_total = [np.nan]

        self.firm_assets_mean = [np.nan]
        self.firm_assets_sd = [np.nan]
        self.firm_assets_total = [np.nan]

        self.firm_deposits_mean = [np.nan]
        self.firm_deposits_sd = [np.nan]
        self.firm_deposits_total = [np.nan]

        self.firm_share_investment_mean = [np.nan]
        self.firm_share_investment_sd = [np.nan]

        self.firm_liabilities_mean = [np.nan]
        self.firm_liabilities_sd = [np.nan]
        self.firm_liabilities_total = [np.nan]

        self.firm_net_worth_mean = [np.nan]
        self.firm_net_worth_sd = [np.nan]
        self.firm_net_worth_total = [np.nan]

        self.firms_final_spot = [np.nan]

        self.bank_assets_mean = [np.nan]
        self.bank_assets_sd = [np.nan]
        self.bank_assets_total = [np.nan]

        self.bank_liabilities_mean = [np.nan]
        self.bank_liabilities_sd = [np.nan]
        self.bank_liabilities_total = [np.nan]

        self.bank_net_worth_mean = [np.nan]
        self.bank_net_worth_sd = [np.nan]
        self.bank_net_worth_total = [np.nan]

        self.firm_dists = []

        # Setup space and firms
        self.product_space = self.setup_product_space(
            self.parameters["number_of_products"],
            self.parameters["product_space_structure"])

        self.prod_space_distances = \
            dict(nx.all_pairs_dijkstra_path_length(self.product_space))

        self.add_complexity_space(
            value_allocation=self.parameters["prod_space_value_allocation"],
            value_distribution=self.parameters["prod_space_value_dist"])

        self.prod_space_vis_position_dict = nx.spring_layout(
            self.product_space, seed=1)

        self.bank_list = self.setup_banks()
        self.firm_list = self.setup_firms()

        self.test_assumptions()

    def test_assumptions(self):
        """Tests model assumptions.
        
        Collects tests for the model assumptions.
        """
        assert len(self.firm_list) == self.n_firms, \
            "Not correct number of firms generated: {} instead of {}".format(
                len(self.firm_list), self.n_firms)
        assert \
            sum([isinstance(i, Firm) for i in self.firm_list]) == self.n_firms,\
            "Not all elements of firm_list are instances of class Firm."
        assert isinstance(self.product_space, nx.classes.graph.Graph), \
            "Product space not network but class {}".format(
                type(self.product_space))
        pass

    def setup_firms(self, allocation_form="small_dc", double_positions=True):
        """Creates firm list and places them on model product space.
        
        Called when the model is initialized.
        The function creates a list of firm instances according to the
        parameters as set out in the parameter file associated with the 
        model.
        The firms are placed on the model-specific product space according
        to a specified 'allocation_form'.
        Whether more than one firm can start on the same spot can also be
        specified via the function arguments.
        
        If an allocation form according to the smallest values of a 
        certain node property (such as centrality) is used, the function 
        proceeds as follows:
            (1) the property is computed for every single node in the 
                product space. This results in a dictionary.
            (2) a  minimum number of possible spots to be given to the 
                firms is specified. This is currently set to 20% of the 
                product space.
            (3) the static method 'extract_small_keys' is used to 
                identify those nodes that have the smallest values for 
                this property. Here, the minimum number of nodes as 
                specified in step (2) is used.
            Finally, for each firm one random position is chosen from the
            resulting list.
        
        Parameters
        ----------
        allocation_form : str, optional
            Specifies how the start positions for the firms on the product 
            space should be chosen, by default "random".
            Currently, the following options are available:
                "random": picks random positions on the product space
                "small_evc": firms are allocated on those nodes that 
                    have the smallest eigenvector centrality (uses 
                    static method `extract_small_keys`)
                "small_dc": firms are allocated on those nodes that 
                    have the smallest degree centrality (uses static 
                    method  `extract_small_keys`)

        double_positions : bool, optional
            Can two or more firms start on the same spot? By default True
        
        Returns
        -------
        firm_list : list
            A list with instances of class firm.Firm
        
        Raises
        ------
        InternalError
            Happens only if there is an internal malfunction and an
            allocation form has been chosen that has not been previously
            defined
        """
        allocation_forms = ("random", "small_evc", "small_dc")
        assert allocation_form in allocation_forms, \
            "Wrong allocation form given, must be one of {} but not {}".format(
                allocation_forms, allocation_form)
        assert isinstance(double_positions, bool), \
            "double_positions must be bool, not {}".format(
                type(double_positions))

        if not double_positions:
            assert self.product_space.number_of_nodes() < self.nb_of_timesteps,\
                "Double pos in product space not allowed, " + \
                "but more firms than products!"

        if allocation_form == "random":
            firm_position = np.random.choice(self.product_space.nodes(),
                                             size=self.n_firms,
                                             replace=double_positions)
        elif allocation_form == "small_evc":
            centralities = nx.eigenvector_centrality(self.product_space)
            # Currently: 20% of the product space
            minimum_nb_positions = \
                int(0.2 * self.product_space.number_of_nodes())
            possible_positions = self.extract_small_keys(
                dict_used=centralities, min_nb_positions=minimum_nb_positions)
            firm_position = np.random.choice(possible_positions.get("indices"),
                                             size=self.n_firms,
                                             replace=double_positions)
        elif allocation_form == "small_dc":
            centralities = nx.degree_centrality(self.product_space)
            # Currently: 20% of the product space
            minimum_nb_positions = \
                int(0.2 * self.product_space.number_of_nodes())
            possible_positions = self.extract_small_keys(
                dict_used=centralities, min_nb_positions=minimum_nb_positions)
            firm_position = np.random.choice(possible_positions.get("indices"),
                                             size=self.n_firms,
                                             replace=double_positions)
        else:
            raise InternalError(
                "Must not happen. If statement for allocation_form is missing!")

        firm_list = [Firm(identifier=i,
                          model_instance=self,
                          init_prod_space_position=firm_position[i - 1])
                     for i in range(1, self.n_firms + 1)]
        return firm_list

    def setup_banks(self):
        """Creates bank list

        Called when the model is initialized.
        The function creates a list of bank instances according to the
        parameters as set out in the parameter file associated with the 
        model.
        """
        bank_list = [Bank(self.financial_regime_parameter,
                          rate_on_deposits=self.rate_on_deposits,
                          rate_on_loans=self.rate_on_loans,
                          payback_rate=self.payback_rate,
                          identifier=i, model_instance=self
                          )
                     for i in range(1, self.n_banks + 1)]

        return bank_list

    def loan_interest_payment(self, firm_id, bank_id):
        """Payment of interest on loans

        Determines amount of payment and respective (paying) firm and (
        receiving) bank.
        Then, calls the method 'payment_firm_to_bank' to execute the
        transaction.

        Parameters
        ----------
        firm_id : int
            The id of the firm that has to pay interests
        bank_id : int
            The id of the bank that receives the interest
        """
        bank_instance = self.bank_list[bank_id - 1]
        interest_payment = bank_instance.get_interest_for_loans(firm_id)
        self.payment_firm_to_bank(firm_id, bank_id, interest_payment)

    def payment_firm_to_bank(self, firm_id, bank_id, payment):
        """Money transaction from firms to banks
        
        The deposits of the firm get reduced, the deposits
        of the bank increased.
        Currently, deposits of the firms must not become negative, 
        defaults are not considered.

        Parameters
        ----------
        firm_id : int
            id of the paying firm
        bank_id : int
            id of the receiving bank
        payment : float
            The payment being processed
        """
        firm_instance = self.firm_list[firm_id-1]
        bank_instance = self.bank_list[bank_id-1]
        # We only have one bank, so choose single account of the firm
        firm_instance.assets["deposits"][bank_id] -= payment
        bank_instance.assets["deposits"] += payment
        bank_instance.liabilities["deposits"][bank_id] += payment

        assert bank_instance.liabilities["deposits"][bank_id] <= 0, \
            "Bank liabilities must be negative!"
        assert firm_instance.assets["deposits"][bank_id] >= 0, \
            "Firm deposits must be positive!"
        assert bank_instance.assets["deposits"] >= 0, \
            "Bank deposits must be positive!"

    def deposit_interest_payment(self, firm_id, bank_id):
        """Payment of interest on deposits

        Determines amount of payment and respective (receiving) firm 
        and (paying) bank. Then, calls the method 'place_deposit' to 
        execute the transaction.

        Parameters
        ----------
        firm_id : int
            The id of the firm that receives interest payments
        bank_id : int
            The id of the bank that has to pay interests
        """
        bank_instance = self.bank_list[bank_id - 1]
        interest_payment = bank_instance.get_interest_on_deposits(firm_id)
        self.place_deposit(firm_id=firm_id, bank_id=bank_id,
                           deposit_value=interest_payment)

    def place_deposit(self, firm_id, bank_id, deposit_value):
        """Money transaction from banks to firms
        
        Adds a deposit for the firm at the bank.
        For the firm this is recorded as an asset, for the bank this is
        recorded as a liability.

        Parameters
        ----------
        firm_id : int
            The id of the firm instance.
        bank_id : int
            The id of the bank instance.
        deposit_value : float
            The amount of the deposit to be created.
        """
        firm_instance = self.firm_list[firm_id-1]
        bank_instance = self.bank_list[bank_id-1]

        if bank_id in firm_instance.assets["deposits"].keys():
            firm_instance.assets["deposits"][bank_id] += deposit_value
        else:
            firm_instance.assets["deposits"][bank_id] = deposit_value

        if firm_id in bank_instance.liabilities["deposits"].keys():
            bank_instance.liabilities["deposits"][firm_id] -= deposit_value
        else:
            bank_instance.liabilities["deposits"][firm_id] = -deposit_value
        if deposit_value < 0:
            pdb.set_trace()
        assert deposit_value >= 0, "deposit_value has wrong sign"

    def place_loan(self, firm_id, bank_id, loan_size):
        """Placement of a loan and money transaction from banks to firms.
        
        Adds a loan.
        For the bank this is recorded as an asset, for the firm it is
        recorded as a liability.
        At the same time, by calling the method 'place_deposit', the 
        loan is transferred to the deposits of the firm, creating a new 
        liability for the bank, and a new asset for the firm.

        Parameters
        ----------
        firm_id : int
            The id of the firm instance.
        bank_id : int
            The id of the bank instance.
        loan_size : float
            The amount of the loan to be created.
        """
        firm_instance = self.firm_list[firm_id-1]
        bank_instance = self.bank_list[bank_id-1]

        # The loan gets added as a liability to the firm:
        if bank_id in firm_instance.liabilities["loans"].keys():
            firm_instance.liabilities["loans"][bank_id] -= loan_size
        else:
            firm_instance.liabilities["loans"][bank_id] = -loan_size

        # The loan gets added as an asset for the bank:
        if firm_id in bank_instance.assets["loans"].keys():
            bank_instance.assets["loans"][firm_id] += loan_size
        else:
            bank_instance.assets["loans"][firm_id] = loan_size
        if loan_size < 0:
            pdb.set_trace()
        assert loan_size >= 0, "Loan has wrong sign."
        # The loan gets transferred on the account of the firm:
        self.place_deposit(firm_id, bank_id, loan_size)

    def repay_loan(self, firm_id, bank_id, payback_amount):
        """Transaction of loan installments from firms to banks

        Current loan installments are paid by firms and received by 
        banks, thereby adjusting their assets and liabilities.

        Parameters
        ----------
        firm_id : int
            The id of the firm instance.
        bank_id : int
            The id of the bank instance.
        payback_amount : float
            The amount of the installment.
        """
        bank_instance = self.bank_list[bank_id-1]
        firm_instance = self.firm_list[firm_id-1]

        # Adjust assets and liabilities of the bank
        bank_instance.receive_loan_payback(firm_id, payback_amount)

        # Adjust assets and liabilities of the firm
        firm_instance.payback_loan(bank_id, payback_amount)

    def setup_product_space(self, n_products, network_structure):
        """Initialize the product space of the model

        Sets up the product space to be used throughout the model run.

        Parameters
        ----------
        n_products : int
            The number of products, i.e. the nb of vertices.
        network_structure : list
            Information about the structure of the network.
            The first element for the list contains a string indicating the
            structure, or a function to be used to the network creation.
            The following elements then specify (if necessary) the
            arguments for this network structure.
            Currently allowed structures:
            ['full']: a full network with each product connected to any 
                other product with the same weight.
            ['ring', n]: a ring network, where each product is connected 
                to n nearest neighbors.
            ['random', p]: an Erdos-Renyi graph with paramter p
            ['BA', m]: A BA graph where m denotes the number of edges to
                       attach from a new node to existing nodes.
            ['powerlaw_cluster', m, p]: A graph with a power law degree
                                        distribution.  m denots the number of
                                        random edges to add for each new node,
                                        p denotes the probability of adding a
                                        triangle after adding a random edge
            ['powerlaw_tree', gamma]: A tree with power law degree distribution,
                                      which has an exponent gamma.

        Raises
        ------
        InputError
            Happens if network structure does not correspond to the currently
            allowed structures.
        """
        assert isinstance(n_products, int), \
            "n_products not int but {}".format(type(n_products))
        if network_structure[0] == "full":
            print(file_marker + "Create full network for the model")
            prod_space = nx.complete_graph(n_products)

        elif network_structure[0] == "ring":
            print(file_marker + "Create ring network for the model")
            prod_space = nx.watts_strogatz_graph(
                n_products, network_structure[1], 0.0)

        elif network_structure[0] == "random":
            print(file_marker + "Create Erdos-Renyi network for the model")
            assert len(network_structure) == 2, \
                "For ER graph exactly two arguments are required, not {}"\
                .format(network_structure)
            prod_space = nx.fast_gnp_random_graph(n_products,
                                                  network_structure[1])

        elif network_structure[0] == "regular":
            print(file_marker + "Create regular network for the model")
            assert len(network_structure) == 2, \
                "For regular graph exactly 2 arguments are required, not {}"\
                .format(network_structure)
            if (network_structure[1] * n_products) % 2 != 0:
                d_parameter = network_structure[1] + 1
                print("Corrected d of regular graph from ",
                      "{} to {}, otherwise n*d would be odd".format(
                          network_structure[1], d_parameter))
            else:
                d_parameter = network_structure[1]
            prod_space = nx.random_regular_graph(d=d_parameter,
                                                 n=n_products)

        elif network_structure[0] == "BA":
            print(file_marker + "Create Barabasi-Albert network for the model")
            # m: Number of edges to attach from a new node to existing nodes
            assert len(network_structure) == 2, \
                "For ER graph exactly two arguments are required, not {}"\
                .format(network_structure)
            assert isinstance(network_structure[1], int), \
                "For BA graph, second argument must be int, not {}".format(
                    type(network_structure[1]))
            prod_space = nx.barabasi_albert_graph(n_products,
                                                  network_structure[1])

        elif network_structure[0] == "powerlaw_cluster":
            print(file_marker + "Create powerlaw cluster network for the model")
            # m: Number of edges to attach from a new node to existing nodes
            # p: prob of transforming triples to triangles
            assert len(network_structure) == 3, \
                "For powerlaw cluster graph 3 arguments are required, " + \
                "not {}".format(network_structure)
            powerlaw_m = network_structure[1]
            assert isinstance(powerlaw_m, int), \
                "For powerlaw cluster graph, 2nd argument must be int, " + \
                "not {}".format(type(powerlaw_m))
            powerlaw_p = network_structure[2]
            assert isinstance(network_structure[2], float) & \
                   (0.0 <= powerlaw_p <= 1.0), \
                   "For powerlaw cluster graph, 3rd arg must be float" + \
                   "betw 0 and 1, not {}".format(
                    type(network_structure[2]))
            # m: the number of random edges to add for each new node
            # p: Probability of adding a triangle after adding a random edge
            prod_space = nx.powerlaw_cluster_graph(n_products,
                                                   powerlaw_m,
                                                   powerlaw_p)

        elif network_structure[0] == "powerlaw_tree":
            print(file_marker + "Create power law tree network for the model")
            # gama: Exponent of the power law.
            assert len(network_structure) == 2, \
                "For powerlaw cluster graph exactly two arguments are " + \
                "required, not {}".format(network_structure)
            prod_space = nx.random_powerlaw_tree(n_products,
                                                 network_structure[1])
        else:
            raise InputError(
                "No correct network structure for the product space given!")

        return prod_space

    def add_complexity_space(self,
                             value_allocation="random",
                             value_distribution={'uniform': [1.0]}):
        """Add information values to all nodes of the product space:
        
        Adds complexity values to the nodes of the product space. 
        First, one has to specify whether complexity values
        are distributed to the nodes randomly or according to a certain 
        node property, such as degree centrality.
        Second, one has to choose the distribution of complexity values, 
        such as uniform, power law or empirical. In the latter case one has
        to provide a list with complexity values. The concrete options are
        described below.
        
        If `value_allocation` is not `random`, then nodes are sorted according
        to the specified property. Then complexity values are sorted as well
        and the node with the highest value for the relevant property gets the
        highest complexity, etc.

        The weights of the edges are set according to the difference of 
        complexity of the products: the larger the different, the greater the
        distance between the products. The products with the largest difference
        in complexity have a distance of 1, for others it is their difference
        divided by the largest difference.
        
        Parameters
        -----------
        value_allocation: str, optional
            Should the complexity values be randomly allocated to the nodes
            of the network, or should the complexity correlate with a 
            property of the node?
            Possible values currently are:
            'random': complexity values are allocated randomly
            'degree_centrality': higher levels of complexity go to nodes with
                a higher degree centrality.
            'eigenvector_centrality': higher levels of complexity go to nodes 
                with a higher eigenvector centrality.
            'closeness_centrality': higher levels of complexity go to nodes 
                with a higher closeness centrality.        
        
        value_distribution: dict, optional
            Gives the kind of distribution according to which complexity values
            for the products will be distributed. The key of the dict should
            indicate the kind of distribution, the value should be a list
            with the parameters. If the key is 'empirical', the value shoud be
            a list with the complexity values to be used.
            The following key-value pairs are allowed:
                'uniform' : [complexity_value]
                'exponential' : [parameters of np.random.exponential]

        Raises
        ------
        InputError
            Happens if value distribution does not correspond to the currently
            possible distributions.

        InternalError
            Happens if a wrong value allocation is given.
        """
        possible_values_value_allocation = ('random', 'degree_centrality',
                                            'eigenvector_centrality',
                                            'closeness_centrality')
        assert value_allocation in possible_values_value_allocation, \
            "Wrong value for value_allocation: {}. Allowed are: {}".format(
                value_allocation, possible_values_value_allocation)

        assert isinstance(value_distribution, dict), \
            "value_distribution should be dict, not {}".format(
                type(value_distribution))

        kind_value_dist = list(value_distribution.keys())
        assert len(kind_value_dist) == 1, \
            "value_distribution must only have one key, but has {}".format(
                value_distribution.keys())

        possible_values_kind_value_dist = \
            ('uniform', 'exponential', 'empirical', "normal")

        if kind_value_dist[0] == 'uniform':
            print(file_marker + \
                  "Values for product complexity created via uniform dist")
            assert isinstance(value_distribution['uniform'], list) & \
                   len(value_distribution['uniform']) == 1, \
                   "Values for uniform dist in value_distribution must be" + \
                   "list with one element, not {}".format(
                    value_distribution['uniform'])
            complexity_values = self.parameters["number_of_products"] \
                * value_distribution['uniform']

        elif kind_value_dist[0] == 'exponential':
            assert isinstance(value_distribution['exponential'], list) & \
                   len(value_distribution['exponential']) == 1, \
                   "Values for exp dist in value_distribution must be" + \
                   "list with one element, not {}".format(
                    value_distribution['exponential'])
            complexity_values = np.random.exponential(
                scale=value_distribution['exponential'],
                size=self.parameters["number_of_products"])

        elif kind_value_dist[0] == 'normal':  # np.random.normal(mu, sigma, 1000)
            assert isinstance(value_distribution['normal'], list), \
                "value_distribution['normal'] must be list not {}".format(
                    type(value_distribution['normal']))
            assert len(value_distribution['normal']) == 2, \
                "Values for normal dist in value_distribution must be list" + \
                "with two elements, not {}".format(
                    value_distribution['normal'])

            complexity_values = np.random.normal(
                loc=value_distribution['normal'][0],
                scale=value_distribution['normal'][1],
                size=self.parameters["number_of_products"])
            complexity_values = [i + abs(min(complexity_values)) + 0.01
                                 for i in complexity_values]

        else:
            error_message = "Wrong specification for complexity_values." + \
                            "Given: {}. But allowed are only: {}".format(
                             kind_value_dist[0],
                             possible_values_kind_value_dist)
            raise InputError(error_message)

        # Now allocate the complexity values to the nodes:
        if value_allocation == "random":
            random_indices = list(range(self.parameters["number_of_products"]))
            np.random.shuffle(random_indices)
            for i in range(self.parameters["number_of_products"]):
                rdn_node_index = random_indices[i]
                self.product_space.nodes[rdn_node_index]["complexity"] = \
                    complexity_values[i]
                self.product_space.nodes[rdn_node_index]["delta"] = \
                    complexity_values[i] * self.delta_coefficient
                self.product_space.nodes[rdn_node_index]["aggregate_output"] = \
                    self.aggregate_output_initial / self.parameters[
                        "number_of_products"]
                self.product_space.nodes[rdn_node_index]["firms"] = 0

        else:
            if value_allocation == "degree_centrality":
                relev_prop = nx.degree_centrality(self.product_space)
            elif value_allocation == "eigenvector_centrality":
                relev_prop = nx.eigenvector_centrality(self.product_space)
            elif value_allocation == "closeness_centrality":
                relev_prop = nx.closeness_centrality(self.product_space)
            else:
                raise InternalError(
                    "No correct value for value_allocation given despite test!")
            sorted_prop = np.asarray(list((
                {k: v for k, v in sorted(relev_prop.items(),
                                         key=lambda item: item[1])}.keys())))
            sc = sorted(complexity_values)
            node_indices = list(range(self.parameters["number_of_products"]))
            if len(set(relev_prop.values())) == 1:
                print("Allocation by {}, yet all nodes have same value!".format(
                    value_allocation))
                np.random.shuffle(node_indices)
            for i in range(self.parameters["number_of_products"]):
                i_node = node_indices[i]
                prop = sorted_prop[i_node]
                self.product_space.nodes[prop]["complexity"] = sc[i]
                self.product_space.nodes[prop]["delta"] = \
                    sc[i] * self.delta_coefficient
                self.product_space.nodes[i_node]["initial_aggregate_output"] = \
                    self.aggregate_output_initial / \
                    self.parameters["number_of_products"]
                self.product_space.nodes[i_node]["aggregate_output"] = 0
                self.product_space.nodes[i_node]["firms"] = 0

        unweighted_weight_list = []
        if kind_value_dist[0] != "uniform":
            for (product1, product2) in self.product_space.edges():
                unweighted_weight = \
                    abs(self.product_space.nodes[product1]["complexity"] -
                        self.product_space.nodes[product2]["complexity"])
                unweighted_weight_list.append(unweighted_weight)
            max_weight = max(unweighted_weight_list)
            for (product1, product2) in self.product_space.edges():
                weighted_weight = \
                    abs(self.product_space.nodes[product1]["complexity"] -
                        self.product_space.nodes[product2]["complexity"]) / \
                    max_weight
                self.product_space.add_weighted_edges_from([(product1, product2,
                                                             weighted_weight)])
        else:
            for (product1, product2) in self.product_space.edges():
                self.product_space.add_weighted_edges_from([(product1, product2,
                                                             1.0)])
        print(file_marker + "Finished allocation of complexity values.",
              "Allocation: {}. Values: {}".format(value_allocation,
                                                  kind_value_dist[0]))

    def visualize_product_space(self, save_fig=True, mark_firms=False,
                                individual_start="", individual_end=""):
        """Visualize the product space of the model

        Creates a visualization of the product space used in this model.

        Parameters
        ----------
        save_fig : bool, optional
            If True the figure will be saved, by default False
        mark_firms : bool, optional
            If true, the position of firms on the product space is marked,
            by default False
        individual_start : str
            Optional prefix for the figure file name.
        individual_end : str
            Optional suffix for the figure file name.

        """
        assert all(map(lambda x: isinstance(x, bool), (save_fig, mark_firms))),\
            "both save_fig and mark_firms should be bool, not {}, {}".format(
                save_fig, mark_firms)
        G = self.product_space
        G_title = "Product space of model " + self.base_name + "_" + \
                  str(self.id) + ": " + str(individual_end)

        fig, ax = plt.subplots(figsize=(8, 8))
        pos = self.prod_space_vis_position_dict
        ax.set_title(G_title)
        ax.axis('off')
        if mark_firms:
            if self.t == -1:
                self.add_firm_positions()
            complexity_list = [50 + (G.nodes[i]["complexity"] * 150)
                               for i in G.nodes]
            firm_list = [G.nodes[i]["firms"] for i in G.nodes]

            nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes,
                                        node_color=firm_list,
                                        node_size=complexity_list,
                                        alpha=0.95, ax=ax)
            nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.15, ax=ax)
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(nc, cax=cbar_ax)
        else:
            nx.draw_networkx_nodes(G, pos, node_size=25, alpha=0.6, ax=ax)
            nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.15, ax=ax)
            plt.tight_layout()

        if save_fig:
            fig_name = str(individual_start) + "product_space_m" + \
                       self.base_name + "_" + str(self.id) + "_" + \
                       str(individual_end) + ".pdf"
            print(file_marker + "Save product space figure in: " + fig_name)
            plt.savefig(fig_name)

    def get_all_visible_products(self, prod_position,
                                 range_of_vision, variant=1):
        """Returns all the products the firm can see
        
        The total number of products the firm can see depends on 
        `range_of_vision`, which is always an integer and specifies the nb
        of products the firm can see.
        The products it can see are the *closest* products. These are not 
        necessarily the immediate neighbors to the current position on the
        product space. Therefore, this function computes the `range_of_vision`
        closest products to the current position.
        
        Note
        -----
        Currently, two variants exist. It is not clear which one is faster.
        Both give the same result.
        
        Parameters
        ----------
        prod_position : int
            the firm's current position on the product space or the product of 
            interest; in both cases vantage point for the search
        range_of_vision : int
            the number of nodes the firm can see
        variant : int
            The computation variant. Variant 1 refers to the matrix of all
            product distances and filters a dictionary to get the closest 
            products. Variant 2 uses `nx.single_source_dijkstra_path_length`
            to compute the nearest products in every time step. So far it is
            unclear which one is faster.

        Returns
        -------
        all_visible_nodes : dict
            a dictionary with all nodes that are visible to the firm and their
            respective distances

        Raises
        ------
        InputError
            Happens if a non-existent variant was chosen.
        """
        if variant == 1:
            range_of_vision += 1
            all_visible_nodes = dict(
                sorted(self.prod_space_distances[prod_position].items(),
                       key=itemgetter(1))[:int(range_of_vision)])
        elif variant == 2:
            all_visible_nodes_ = nx.single_source_dijkstra_path_length(
                G=self.product_space,
                source=prod_position,
                cutoff=range_of_vision
            )
            range_of_vision += 1
            all_visible_nodes = dict(
                sorted(all_visible_nodes_.items(),
                       key=itemgetter(1))[:range_of_vision])
        else:
            raise InputError(
                "False way to compute visible products chosen: {}".format(
                    variant))

        return all_visible_nodes

    def run(self):
        """Runs the model

        1. Updates firm positions on the product space
        2. Computes the prices for all products on the product space
        3. Compute output of each firm
            After computing the potential output, for every product it is 
            checked whether the aggregated output exceeds the nominal demand.
            If this is the case, each firm gets a share of the nominal demand 
            according to its share in the potential total output of all firms.
        4. Firms invest and update capital stock using `Firm.update_firm()`
            They (1) compute inventories, (2) choose their target node, 
            (3) compute actual profits, (4) fix desired investment and apply 
            for credit, (5) conduct R&D, (6) actually move on the product space
            and (7) update their properties (account, capital stock, etc.)
        5. Update the bank accounts
            Just calls the method `update_bank` from the bank instance
        6. Updates the output on the product space
        7. Computes and saves the relevant summary statistics
        """
        # get aggregate output
        print("------Called run function-------------------------------------")
        record_distributions = (1, self.nb_of_timesteps - 1)
        for i in range(self.nb_of_timesteps):
            print(file_marker + "Model {}: t {}/{}".format(
                self.id, i, self.nb_of_timesteps), end="\r")
            self.t = i
            # 1. Update firm positions on the product space
            self.add_firm_positions()

            # 2. Compute the prices for all products on the product space
            self.add_prices()  # If output=0, then delta is used as price

            # 3.1 Reset production statistics on product space:
            old_output_dict = nx.get_node_attributes(
                self.product_space, "aggregate_output")

            excess_nodes = \
                {k: v for k, v in old_output_dict.items() \
                 if v > self.nominal_demand}
            assert len(excess_nodes.keys()) == 0, \
                "For some nodes production > nominal demand: {}".format(
                    excess_nodes)

            nx.set_node_attributes(self.product_space, 0.0, "aggregate_output")

            new_output_dict = {}
            capital_stock_dict = {}
            all_deposits = 0
            all_loans = 0
            individual_output_dict = {}
            total_output_dict = dict(  # Starts with value 0 for each node
                zip(self.product_space.nodes(),
                    self.parameters["number_of_products"] * [0]))

            # 3.2 Compute output for each firm, save it,
            # add it to total_output_dict
            for firm in self.firm_list:
                price = self.product_space.nodes[ \
                    firm.product_space_position]["price"]
                capital_stock = firm.capital_stock[-1]
                output = Firm.output_and_capability_costs(
                    firm, capital_stock, price)
                if math.isnan(output) == True:
                    pdb.set_trace()
                total_output_dict[firm.product_space_position] += output
                individual_output_dict[firm] = output

            # 3.3 If total output for a product exceeds demand,
            # update actual outputs
            nodes_with_excess_demands = \
                {k: v for k, v in total_output_dict.items() \
                 if v > self.nominal_demand}

            for k in nodes_with_excess_demands.keys():
                for firm in self.firm_list:
                    price = self.product_space.nodes[ \
                        firm.product_space_position]["price"]
                    if firm.product_space_position == k:
                        possible_output = \
                            firm.output_and_capability_costs(
                                firm.capital_stock[-1], price)
                        adapted_output = \
                            possible_output / total_output_dict[k] \
                            * self.nominal_demand
                        if math.isnan(adapted_output):
                            pdb.set_trace()
                        individual_output_dict[firm] = \
                            self.round_half_down(adapted_output, 2)

            # 4. Location decision, investment of the firms and
            # update capital stock
            for firm in self.firm_list:
                # 4.0 Save preparatory values
                output = individual_output_dict.get(firm)
                total_output_this_period = \
                    total_output_dict.get(firm.product_space_position)
                price = self.product_space.nodes[
                    firm.product_space_position]["price"]
                total_output_prev_t = \
                    old_output_dict[firm.product_space_position]

                # 4.1 Write firm output into new_output_dict
                if firm.product_space_position in new_output_dict.keys():
                    new_output_dict[firm.product_space_position] = \
                        int(new_output_dict[firm.product_space_position] +
                            output)
                else:
                    new_output_dict[firm.product_space_position] = output

                assert new_output_dict[firm.product_space_position] <= \
                       self.nominal_demand, \
                       "new_output_dict value {} exceeds maximum {}".format(
                        new_output_dict[firm.product_space_position],
                        self.nominal_demand)

                # 4.2 Firms choose new spot on the product space, compute actual
                #   profits, make investment decisions demand money from banks 
                #   and update capital account
                capital_stock, account_firm = firm.update_firm(
                    price=price,
                    financial_regime_parameter=self.financial_regime_parameter,
                    output=output,
                    total_output_this_period=total_output_this_period)
                if account_firm > 0:
                    all_deposits += account_firm
                else:
                    all_loans += account_firm
                capital_stock_dict[firm] = capital_stock

            # 5. Update bank accounts
            account_bank = self.bank_list[0].update_bank(all_deposits,
                                                         all_loans)
            assert account_bank == all_loans + all_deposits, \
                "ALERT: account_bank={}, all_loans + all_deposits={}+{}={}!"\
                    .format(account_bank,
                            all_loans,
                            all_deposits,
                            all_loans + all_deposits)

            # 6. Update the output for each product *on the product space*
            excess_nodes = \
                {k: v for k, v in new_output_dict.items()
                 if v > self.nominal_demand}
            assert len(excess_nodes.keys()) == 0, \
                "For some nodes production > nominal demand ({}): {}. " \
                "Check also nodes_with_excess_demands: {}".format(
                    self.nominal_demand,
                    excess_nodes,
                    nodes_with_excess_demands)

            nx.set_node_attributes(self.product_space, new_output_dict,
                                   "aggregate_output")

            # 7. Compute relevant statistics 
            # Get the complexities of those products actually produced:
            comp_prod_prdcts = \
                [self.product_space.nodes[i]["complexity"] for i in \
                 self.product_space.nodes if \
                 self.product_space.nodes[i]["aggregate_output"] > 0.0]

            assert len(comp_prod_prdcts) > 0, "NOTHING GETS PRODUCED!!!"
            assert np.isnan(comp_prod_prdcts).any() == False, \
                "Nan values in comp_prod_prdcts ({})".format(comp_prod_prdcts)

            # Get the prices of those products actually produced:
            price_prod_prdcts = \
                [self.product_space.nodes[i]["price"] for i in \
                 self.product_space.nodes if \
                 self.product_space.nodes[i]["aggregate_output"] > 0.0]
            assert len(price_prod_prdcts) > 0, "NOTHING HAS A PRICE!!!"
            assert np.isnan(price_prod_prdcts).any() == False, \
                "Nan values in price_prod_prdcts ({})".format(price_prod_prdcts)

            self.record_state_variables(
                time=i,
                new_output_dict=new_output_dict,
                comp_prod_prdcts=comp_prod_prdcts,
                capital_stock_dict=capital_stock_dict,
                price_prod_prdcts=price_prod_prdcts,
                t_record_distributions=record_distributions)

    def record_state_variables(self, time, new_output_dict,
                               comp_prod_prdcts, capital_stock_dict,
                               price_prod_prdcts,
                               t_record_distributions):
        """Record state variables of the model
        
        Used to save all the state variables of interest.
        
        Parameters
        ----------
        time : int
            Current time step
        new_output_dict : dict
            Dict with keys being produced on the product space and values
            the total value produced of these products
        comp_prod_prdcts : list
            List of complexity values of all products that have output > 0
        capital_stock_dict : dict
            Dict with firms as keys and their capital stock as value
        price_prod_prdcts : list
            List with the prices of all products that have output > 0
        t_record_distributions : list
            List including those time steps in which distributions should be
            saved using the method `self.record_distributions()`
        """
        self.agg_output_all_prdcts.append(
            sum([v for v in new_output_dict.values()]))

        self.products_produced_share.append(
            len(comp_prod_prdcts) / self.product_space.number_of_nodes())
        self.agg_capital_stock.append(sum([w for w in
                                           capital_stock_dict.values()]))

        self.products_produced_avg_comp.append(np.mean(comp_prod_prdcts))
        self.prize_avg.append(np.mean(price_prod_prdcts))
        self.prize_sd.append(np.std(price_prod_prdcts))
        try:
            self.prize_compl_cor.append(
                pearsonr(comp_prod_prdcts, price_prod_prdcts)[0])
        except ValueError:
            self.prize_compl_cor.append(np.nan)

        firm_innovation_capabilities = []
        firm_spillover_capabilities = []
        firm_rd_investment = []
        firm_ac_investment = []
        firm_k_investment = []
        firm_profits = []
        firm_gross_revenues = []
        firm_accounts = []
        firm_assets = []
        firm_deposits = []
        firm_liabilities = []
        firm_net_worth = []
        bank_assets = []
        bank_liabilities = []
        bank_net_worth = []
        firms_final_position = []
        firm_share_investment = []

        for firm in self.firm_list:
            firm_innovation_capabilities.append(
                firm.innovation_capabilities[-1])
            firm_spillover_capabilities.append(firm.spillover_capabilities[-1])
            firm_rd_investment.append(firm.rd_investment)
            firm_ac_investment.append(firm.ac_investment)
            firm_k_investment.append(firm.k_investment)
            firm_profits.append(firm.profit[-1])
            firm_gross_revenues.append(firm.gross_revenue[-1])
            firm_accounts.append(firm.account_firm[-1])
            firm_assets.append(firm.firm_assets)
            firm_liabilities.append(firm.firm_liabilities)
            firm_net_worth.append(firm.firm_net_worth)
            firm_deposits.append(firm.firm_deposits)
            firms_final_position.append(firm.final_spot)
            firm_share_investment.append(firm.share_deposits_invested)
        for bank in self.bank_list:
            bank_balance_sheet = bank.get_balance_sheet()
            bank_assets.append(bank_balance_sheet["assets"])
            bank_liabilities.append(bank_balance_sheet["liabilities"])
            bank_net_worth.append(bank_balance_sheet["net_worth"])

        self.firm_innovation_capabilities_mean.append(
            np.mean(firm_innovation_capabilities))
        self.firm_innovation_capabilities_sd.append(
            np.std(firm_innovation_capabilities))
        self.firm_innovation_capabilities_total.append(
            sum(firm_innovation_capabilities))

        self.firm_spillover_capabilities_mean.append(
            np.mean(firm_spillover_capabilities))
        self.firm_spillover_capabilities_sd.append(
            np.std(firm_spillover_capabilities))
        self.firm_spillover_capabilities_total.append(
            sum(firm_spillover_capabilities))

        self.firm_rd_investment_mean.append(np.mean(firm_rd_investment))
        self.firm_rd_investment_sd.append(np.std(firm_rd_investment))
        self.firm_rd_investment_total.append(sum(firm_rd_investment))

        self.firm_k_investment_mean.append(np.mean(firm_k_investment))
        self.firm_k_investment_sd.append(np.std(firm_k_investment))
        self.firm_k_investment_total.append(sum(firm_k_investment))

        self.firm_ac_investment_mean.append(
            np.mean(firm_ac_investment))
        self.firm_ac_investment_sd.append(
            np.std(firm_ac_investment))
        self.firm_ac_investment_total.append(
            sum(firm_ac_investment))

        self.firm_profits_mean.append(np.mean(firm_profits))
        self.firm_profits_sd.append(np.std(firm_profits))
        self.firm_profits_total.append(sum(firm_profits))

        self.firm_gross_revenues_mean.append(np.mean(firm_gross_revenues))
        self.firm_gross_revenues_sd.append(np.std(firm_gross_revenues))
        self.firm_gross_revenues_total.append(sum(firm_gross_revenues))

        self.firm_accounts_mean.append(np.mean(firm_accounts))
        self.firm_accounts_sd.append(np.std(firm_accounts))
        self.firm_accounts_total.append(sum(firm_accounts))

        self.firm_assets_mean.append(np.mean(firm_assets))
        self.firm_assets_sd.append(np.std(firm_assets))
        self.firm_assets_total.append(sum(firm_assets))

        self.firm_deposits_mean.append(np.mean(firm_deposits))
        self.firm_deposits_sd.append(np.std(firm_deposits))
        self.firm_deposits_total.append(sum(firm_deposits))

        self.firm_liabilities_mean.append(np.mean(firm_liabilities))
        self.firm_liabilities_sd.append(np.std(firm_liabilities))
        self.firm_liabilities_total.append(sum(firm_liabilities))

        self.firm_net_worth_mean.append(np.mean(firm_net_worth))
        self.firm_net_worth_sd.append(np.std(firm_net_worth))
        self.firm_net_worth_total.append(sum(firm_net_worth))

        self.firm_share_investment_mean.append(np.nanmean(firm_share_investment))
        self.firm_share_investment_sd.append(np.nanstd(firm_share_investment))

        self.firms_final_spot.append(sum(firms_final_position) / self.n_firms)

        self.bank_assets_mean.append(np.mean(bank_assets))
        self.bank_assets_sd.append(np.std(bank_assets))
        self.bank_assets_total.append(sum(bank_assets))

        self.bank_liabilities_mean.append(np.mean(bank_liabilities))
        self.bank_liabilities_sd.append(np.std(bank_liabilities))
        self.bank_liabilities_total.append(sum(bank_liabilities))

        self.bank_net_worth_mean.append(np.mean(bank_net_worth))
        self.bank_net_worth_sd.append(np.std(bank_net_worth))
        self.bank_net_worth_total.append(sum(bank_net_worth))

        self.share_firms_waiting.append(sum(
            [f.wait_for_move for f in self.firm_list]
            ) / len(self.firm_list))
        assert 0 <= self.share_firms_waiting[-1] <= 1, \
            "Share of waiting firms: {}".format(self.share_firms_waiting[-1])

        self.share_firms_remained.append(sum(
            [f.remained_same_node for f in self.firm_list]
            ) / len(self.firm_list))
        assert 0 <= self.share_firms_remained[-1] <= 1, \
            "Share of remained firms: {}".format(self.share_firms_remained[-1])

        if time in t_record_distributions:
            self.firm_dists.append(self.record_distributions(
                kind="firms", time=time))

    def record_distributions(self, kind, time):
        """Saves information about relevant distributions of the model

        Parameters
        ----------
        kind : str
            The distribution of what should be visualized; currently 
            only 'firms' is allowed.
        time : int
            Time step in which distributions are saved
        """
        if kind == "firms":
            t_len = len(self.firm_list)
            firm_dist_dict = {
                "firm_innovation_capabilities": \
                    [f.innovation_capabilities[-1] for f in self.firm_list],
                "firm_spillover_capabilities": \
                    [f.spillover_capabilities[-1] for f in self.firm_list],
                "firm_capital_stock": [f.capital_stock[-1]
                                       for f in self.firm_list],
                "firm_profits": [f.profit for f in self.firm_list],
                "id": [self.id] * t_len,
                "t": [time] * t_len
            }
            return_frame = pd.DataFrame.from_dict(firm_dist_dict)
        else:
            raise InputError("Wrong distribution kind given!")
        return return_frame

    def add_firm_positions(self):
        """Update new firm positions

        Update information on where each firm is positioned in the product
        space.
        Called in method 'run' (i.e. once in every time step).
        """
        list_of_positions = [f.product_space_position for f in self.firm_list]
        firm_pos_counter = Counter(list_of_positions)
        capital_dict = {prod: 0 for prod in range(len(self.product_space))}
        for f in self.firm_list:
            capital_dict[f.product_space_position] += f.capital_stock[-1]

        assert sum(firm_pos_counter.values()) == self.n_firms, \
            "Nb of firms on prod space ({}) not equal to firm nb ({})".format(
                sum(firm_pos_counter.values()), self.n_firms)
        nx.set_node_attributes(self.product_space, 0, "firms")
        nx.set_node_attributes(self.product_space, firm_pos_counter, "firms")
        nx.set_node_attributes(self.product_space, capital_dict, "capital")

    def add_prices(self):
        """Update current prices

        Update information on the current price of each product in the
        product space.
        Called once per time step.
        If the product is currently not produced and a firm needs to assume a
        price to compute expected profits, then the monopoly price 'delta' will
        be used.
        """
        for i in range(self.parameters["number_of_products"]):
            self.product_space.nodes[i]["price"] = \
                self.product_space.nodes[i]["delta"] / max(
                    self.product_space.nodes[i]["firms"], 1)

    def return_results(self):
        """
        Returns a pd.DataFrame with the results, such that they can be
        concatenated with results from other simulation runs.
        """
        print(file_marker + "called return_results()")
        ts_len = self.nb_of_timesteps + 1

        result_dict = {
            "id": [self.id] * (self.nb_of_timesteps + 1),
            "t": np.arange(0, self.nb_of_timesteps + 1, 1),
            "n_firms": [self.n_firms] * ts_len,
            "n_workers": [self.n_workers] * ts_len,
            "n_products": [self.parameters["number_of_products"]] * ts_len,
            "prod_space_struc": \
                ['_'.join(str(i) for i in self.parameters[ \
                    "product_space_structure"])] * ts_len,
            "prod_sp_val_alct": [self.parameters[ \
                                     "prod_space_value_allocation"]] * ts_len,
            "prod_sp_val_dist": [str(list(
                self.parameters["prod_space_value_dist"].keys())[0]) + "_" + \
                                 '_'.join(str(e) for e in self.parameters[
                                     "prod_space_value_dist"].values())] * ts_len,
            "delta_coefficient": [self.delta_coefficient] * ts_len,
            "max_info": [self.parameters["max_info"]] * ts_len,
            "initial_capital_stock": [self.initial_capital_stock] * ts_len,
            "initial_capabilities": [self.initial_capabilities] * ts_len,
            "nominal_demand": [self.nominal_demand] * ts_len,
            "fin_regime": [self.financial_regime_parameter] * ts_len,
            "cap_prod": [self.capital_productivity] * ts_len,
            "max_rd_success": [self.max_rd_success] * ts_len,
            "depreciation_rate": [self.depreciation_rate] * ts_len,
            "output_all_prdcts": self.agg_output_all_prdcts,
            "share_produced_products": self.products_produced_share,
            "comp_produced_prod_mean": self.products_produced_avg_comp,
            "price_produced_prod_mean": self.prize_avg,
            "price_produced_prod_std": self.prize_sd,
            "correlation_price_complexity": self.prize_compl_cor,
            "agg_capital_stock": self.agg_capital_stock,
            "firm_inno_cap_mean": self.firm_innovation_capabilities_mean,
            "firm_inno_cap_sd": self.firm_innovation_capabilities_sd,
            "firm_inno_cap_total": self.firm_innovation_capabilities_total,
            "firm_spillover_cap_mean": self.firm_spillover_capabilities_mean,
            "firm_spillover_cap_sd": self.firm_spillover_capabilities_sd,
            "firm_spillover_cap_total": self.firm_spillover_capabilities_total,
            "firm_rd_invest_mean": self.firm_rd_investment_mean,
            "firm_rd_invest_sd": self.firm_rd_investment_sd,
            "firm_rd_invest_total": self.firm_rd_investment_total,
            "firm_ac_invest_mean": self.firm_ac_investment_mean,
            "firm_ac_invest_sd": self.firm_ac_investment_sd,
            "firm_ac_invest_total": self.firm_ac_investment_total,
            "firm_k_invest_mean": self.firm_k_investment_mean,
            "firm_k_invest_sd": self.firm_k_investment_sd,
            "firm_k_invest_total": self.firm_k_investment_total,
            "firm_profits_mean": self.firm_profits_mean,
            "firm_profits_sd": self.firm_profits_sd,
            "firm_profits_total": self.firm_profits_total,
            "firm_gross_revenues_mean": self.firm_gross_revenues_mean,
            "firm_gross_revenues_sd": self.firm_gross_revenues_sd,
            "firm_gross_revenues_total": self.firm_gross_revenues_total,
            "firm_accounts_mean": self.firm_accounts_mean,
            "firm_accounts_sd": self.firm_accounts_sd,
            "firm_accounts_total": self.firm_accounts_total,
            "firm_share_waiting": self.share_firms_waiting,
            "firm_share_remained": self.share_firms_remained,
            "firm_deposits_mean": self.firm_deposits_mean,
            "firm_deposits_sd": self.firm_deposits_sd,
            "firm_assets_mean": self.firm_assets_mean,
            "firm_assets_sd": self.firm_assets_sd,
            "firm_assets_total": self.firm_assets_total,
            "firm_liabilities_mean": self.firm_liabilities_mean,
            "firm_liabilities_sd": self.firm_liabilities_sd,
            "firm_liabilities_total": self.firm_liabilities_total,
            "firm_net_worth_mean": self.firm_net_worth_mean,
            "firm_net_worth_sd": self.firm_net_worth_sd,
            "firm_net_worth_total": self.firm_net_worth_total,
            "bank_assets_mean": self.bank_assets_mean,
            "bank_assets_sd": self.bank_assets_sd,
            "bank_assets_total": self.bank_assets_total,
            "bank_liabilities_mean": self.bank_liabilities_mean,
            "bank_liabilities_sd": self.bank_liabilities_sd,
            "bank_liabilities_total": self.bank_liabilities_total,
            "bank_net_worth_mean": self.bank_net_worth_mean,
            "bank_net_worth_sd": self.bank_net_worth_sd,
            "bank_net_worth_total": self.bank_net_worth_total,
            "firms_final_spot": self.firms_final_spot,
            "firm_share_investment_mean": self.firm_share_investment_mean,
            "firm_share_investment_sd": self.firm_share_investment_sd
        }
        result_frame = pd.DataFrame.from_dict(result_dict)

        results_dists = pd.concat(self.firm_dists, ignore_index=True)
        return_dict = {"dynamics": result_frame,
                       "distributions": results_dists}
        return return_dict

    @staticmethod
    def extract_small_keys(dict_used, min_nb_positions):
        """Extract keys with the smallest corresponding values
        
        Takes a dictionary and an integer n. Returns those keys of the 
        dictionary of which the values are among the n smallest. If there
        are many keys with the same values, more than n keys are returned,
        i.e. those of which the value is smaller or equal to the n-smallest
        value.
        
        Parameters
        ----------
        dict_used : dict
            The dictionary of which the keys are to be extracted
        min_nb_positions : int
            The minimum number of keys to be extracted
        
        Returns
        -------
        dict
            A dictionary with two entries: 'indices', a list with the 
            keys of which the values are smaller or equal to the n smallest
            values; 'threshold_value', a float specifying the threshold value
            used
        """
        assert isinstance(dict_used, dict), \
            "dict_used must be dict, not {}".format(type(dict_used))
        assert isinstance(min_nb_positions, int), \
            "min_nb_positions must be int, not {}".format(type(min_nb_positions))
        sorted_dict = {k: v for k, v in sorted(dict_used.items(),
                                               key=lambda item: item[1])}
        index_value_min_value_pos = list(sorted_dict.keys())[min_nb_positions]
        min_value = sorted_dict[index_value_min_value_pos]
        possible_keys = list(
            {k: v for k, v in sorted_dict.items() if v <= min_value})
        return {"indices": possible_keys, "threshold_value": min_value}

    @staticmethod
    def round_half_down(n, decimals=0):
        """Rounds down a number
                
        Parameters
        ----------
        n : float
            The value to be rounded
        decimals : int, optional
            The nb of decimal places to be rounded to, by default 0
        
        Returns
        -------
        float
            The rounded value
        """
        multiplier = 10 ** decimals
        return math.ceil(n * multiplier - 0.5) / multiplier


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class InternalError(Error):
    """Exception raised internal code inconsistencies.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
