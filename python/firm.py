import numpy as np
import os
from scipy.stats import bernoulli
import networkx as nx
import pdb
import math

file_marker = "[" + str(os.path.basename(__file__)) + "]: "


class Firm:
    """A firm

    Firms are product-seeking agents.
    (1) They produce a product that corresponds to their current 
    position on the product space
    (2) They choose a target-node: within their range of vision firms will
    choose the product that offers highest expected profits as a target
    (3) In order to be able to change their production to the target,
    firms invest into capability-enhancing measures 
    (R&D, absorptive capacities)
    (4) Firms invest into their capital stock according to their 
    production plan
    (5) If capability measures were successful, they change their 
    production.
    """

    def __init__(self, identifier, model_instance,
                 init_prod_space_position, print_progress=False):
        """Initializes a firm instance

        Sets relevant parameters as firm properties.

        Parameters
        ----------
        identifier : int
            A numerical identifier for the firm.
        model_instance : model.Model
            The instance of the associated model
        init_prod_space_position : int
            The initial position on the product space of the model
        print_progress : bool
            Should ongoing information be printed during model execution?
        """
        self.id = identifier
        self.model_instance = model_instance
        self.product_space_position = init_prod_space_position
        self.path = [init_prod_space_position]

        # Initialization values
        self.age = 0
        self.employees = set()
        self.savings = 0
        self.rd_stock = 0
        self.profit = []
        self.gross_revenue = []
        self.account_firm = [0]

        # Currently firms stick to their one bank
        self.firm_bank_id = self.model_instance.bank_list[
            np.random.randint(
                low=0, high=len(self.model_instance.bank_list), size=1)[0]
        ].identifier

        self.assets = {
            "capital_stock": model_instance.parameters["initial_capital_stock"],
            "deposits": {}
        }
        self.liabilities = {
            "loans": {}
        }
        self.firm_net_worth = 0.0
        self.firm_assets = 0.0
        self.firm_liabilities = 0.0
        self.firm_deposits = 0.0

        self.output = []
        self.inventories = [0]
        self.rd_investment = np.nan
        self.ac_investment = np.nan
        self.k_investment = np.nan
        self.wait_for_move = False
        self.remained_same_node = False
        self.available_steps = 1
        self.ultimate_target_node = False

        # Parameters
        params = self.model_instance.parameters
        self.range_of_vision = 1

        self.max_rd_success = params["max_rd_success"]
        self.capital_stock = [params["initial_capital_stock"]]
        self.innovation_capabilities = [params["initial_capabilities_equal"]]
        self.spillover_capabilities = [params["initial_capabilities_equal"]]
        self.depreciation_rate = params["depreciation_rate"]
        self.capital_productivity = params["capital_productivity"]
        self.cost_share = self.capital_productivity / 2
        self.target_prob = params["target_prob"]
        self.p_innovation = params["p_innovation"]
        self.p_spillovers = params["p_spillovers"]
        self.p_information = params["p_information"]
        self.rate_on_deposits = params["rate_on_deposits"]
        self.rate_on_loans = params["rate_on_loans"]
        self.payback_rate = params["payback_rate"]
        self.nominal_demand = params["nominal_demand"]
        self.effect_rd_success = params["effect_rd_success"]
        self.effect_spillover_success = params["effect_spillover_success"]

        self.print_progress = print_progress
        self.firms_in_market = []
        self.price = []
        self.total_output_prod_produced = []
        self.demand_rd = []
        self.total_deposits = [0]
        self.share_deposits_invested = [0]
        self.max_info = params["max_info"]
        self.cant_go_back = params["cant_go_back"]
        self.final_spot = False

    def __repr__(self):
        return 'Firm {} (age: {}, account: {}) at {}'.format(
            self.id, self.age, self.account_firm[-1], hex(id(self))
        )

    def check_final_spot(self):
        """Checks whether the firm has reached its local optimum
        
        Should serve as an indicator of whether the firm can be expected
        to move again in the immediate future. To this end, we check
        whether there is any project on the entire product space for 
        which the firm would currently expect a higher payoff than for
        its current market. In this computation of expected profits
        we abstract from costs and interest payments since they were 
        independent from the product the firm is producing.
        """
        visible_products = self.model_instance.product_space.nodes
        expected_profit_dict = {}
        for prod in visible_products:
            price_p = self.model_instance.product_space.nodes[prod][
                "price"]
            output_p = min(
                self.model_instance.product_space.nodes[prod]["capital"] *
                self.capital_productivity,
                self.nominal_demand)

            if output_p == 0:
                market_share = 1
            else:
                market_share = \
                    self.capital_stock[-1] / \
                    (self.model_instance.product_space.nodes[prod]["capital"] +
                     self.capital_stock[-1])
            output = max(
                min(self.nominal_demand * market_share,
                    self.capital_stock[-1] * self.capital_productivity),
                0)

            expected_profit_dict[prod] = price_p * output

        personal_optimum = max(expected_profit_dict,
                               key=expected_profit_dict.get)

        self.final_spot = self.product_space_position == personal_optimum

    def get_balance_sheet(self):
        """Returns information about the firm's balance sheet
        
        Note: during execution, `self.total_deposits` get updated.
        
        Returns
        -------
        balance_sheet: dict
            Dict with keys `net_worth`, `deposits`, `assets`, and
            `liabilities`.
        """
        total_deposits = self.assets["deposits"][self.firm_bank_id]
        total_assets = self.assets["capital_stock"] + total_deposits
        total_liabilities = sum(
            [sum(list(self.liabilities[i].values())) \
             for i in self.liabilities.keys()])

        assert total_assets >= 0, "Total assets of firm are negative"
        assert total_liabilities <= 0, \
            "Total liabilities of firm are positive"

        net_worth = total_assets + total_liabilities
        balance_sheet = {
            "net_worth": net_worth,
            "deposits": total_deposits,
            "assets": total_assets,
            "liabilities": total_liabilities
        }
        self.total_deposits.append(total_deposits)
        return balance_sheet

    def payback_loan(self, bank_id, payback_amount):
        """Pays back a loan installment
        
        Reduces the open loans by the amount of `payback_amount`. Gets
        called from within the model, so the accompanying adjustment
        of the bank balance sheet is automatically accounted for there.

        Parameters
        ----------
        bank_id : int
            Identifier that corresponds to the bank that receives the
            installment payment.

        payback_amount : float
            The amount the firm pays back.
        """
        self.liabilities["loans"][bank_id] -= payback_amount
        assert self.liabilities["loans"][bank_id] <= 0, \
            "Loans cannot be positive!"

    def pay_receive_interest(self):
        """Pays interest on loans and receives interest for deposits

        Called once per time step per firm in method 'update_firm'
        """
        for b in self.assets["deposits"].keys():
            self.model_instance.deposit_interest_payment(
                firm_id=self.id, bank_id=b)
        for b in self.liabilities["loans"].keys():
            self.model_instance.loan_interest_payment(
                firm_id=self.id, bank_id=b
            )
        if self.assets["deposits"][self.firm_bank_id] < 0:
            pdb.set_trace()
        assert self.assets["deposits"][self.firm_bank_id] >= 0, \
            "Deposits after investments should be >= 0, not {}".format(
                self.assets["deposits"])

    def update_firm(self, price, financial_regime_parameter, output,
                    total_output_this_period):
        """Update the single firm
        
        Goes through the following steps:
        1. Compute the potential output
            Firms aim to produce no more than the theoretical maximum of 
            the nominal demand for a given product.
        2. Compute inventories as the difference between potential and 
            actual output
        3. Choose the target node using the function `choose_target_node`.
            To this end the firm compares expected profits among all 
            products in its reach and chooses the node that promises 
            highest profits as the target node.
        4. Compute actual profits, price_cost_ratio and revenue
        5. Desired investment is fixed
        6. Credit is requested from the bank
        7. Compute financial constraint and fix investment
        8. Conduct capability enhancing R&D
        9. (If possible) go to new spot on the product space
        10. Update capital stock and account
        11. Update firm properties
        
        Parameters
        ----------
        price : float
            Market price of the product.
        financial_regime_parameter : float
            Some parameter > 1 that is used to compute a firm's financial
            constraint. The higher a firm's profit rate, the greater it's
            financial constraint, meaning that it can invest more.
        output : float
            The actual output of the firm in the current time step
            (this has already been updated according to its market share)
        total_output_this_period : float
            The total output in the previous time step for the product 
            the firm is currently producing.

        Returns
        -------
        updated_capital_stock, updated_account_firm : float
            Used in the model class to compute aggregate output
        """
        self.output.append(output)
        self.total_output_prod_produced.append(total_output_this_period)

        # 1. Compute the potential output
        potential_output = self.output_and_capability_costs(
            previous_capital_stock=self.capital_stock[-1],
            price=price)

        # 2. Compute the inventories
        old_inventory = max(self.inventories[-1] - output, 0)
        new_inventory = max((potential_output - output), 0)
        self.inventories.append(old_inventory + new_inventory)

        # 3. Choose the target node
        neighbors = list(self.model_instance.product_space.neighbors(
            self.product_space_position))

        target_node, same_node, target_complexity, target_firms = \
            self.choose_target_node(
                price=price,
                output=output,
            )

        # 4. Compute gross revenues and related values
        price_cost_ratio, profit, payback, gross_revenue, interest_payments = \
            self.get_profit(price=price, output=output)
        if gross_revenue < 0:
            print("gross_revenue<0!!! Aborting!!!")
            pdb.set_trace()

        # Firms put their gross revenue on their bank account,
        # when revenue is negative, they take up a loan to fix it:
        if gross_revenue < 0:
            print("Warning! gross_revenue < 0")
            self.model_instance.place_loan(
                firm_id=self.id, bank_id=self.firm_bank_id,
                loan_size=abs(gross_revenue))
        else:
            self.model_instance.place_deposit(
                firm_id=self.id, bank_id=self.firm_bank_id,
                deposit_value=gross_revenue)

        # Firms pay back loans
        if payback > 0:
            self.model_instance.repay_loan(self.id, self.firm_bank_id, payback)

        # Firms pay interest rates on their loans and receive interest on
        # their deposits
        self.pay_receive_interest()

        # Conclude by saving the state variables
        self.gross_revenue.append(gross_revenue)
        self.profit.append(profit)
        self.price.append(price)
        self.firms_in_market.append(self.model_instance.product_space.nodes[
            self.product_space_position].get(
            "firms"))

        # 5. Compute desired investment
        demand_investment, demand_credit, demand_info, \
        demand_spillovers, demand_k, demand_rd = \
            self.get_demand_for_investment(
                output=output,
                deposits=self.total_deposits[-1],
                target_complexity=target_complexity,
                target_firms=target_firms
            )

        # 6. Request credit from the bank
        # If loan is granted, it's already booked on firm's deposits
        credit_granted = self.model_instance.bank_list[0].credibility(
            firm_id=self.id,
            profit=self.profit[-1],
            capital_stock=self.capital_stock[-1],
            credit_demand=demand_credit,
            deposits=self.total_deposits[-1])

        # 7. Compute financial constraint and fix investment
        financial_constraint = self.assets["deposits"][self.firm_bank_id]
        rd_investment, info_investment, spillover_investment, \
        k_investment, remaining_profit = \
            self.get_investment(
                financial_constraint=financial_constraint,
                demand_k=demand_k,
                demand_rd=demand_rd,
                demand_info=demand_info,
                demand_spillovers=demand_spillovers
            )

        # 8. Conduct capability enhancing R&D
        # Reduce firm deposits by amount of investments:
        total_investments = math.floor(rd_investment + info_investment +
                                       spillover_investment + k_investment)
        if total_investments > self.assets["deposits"][self.firm_bank_id]:
            pdb.set_trace()
        if self.assets["deposits"][self.firm_bank_id] == 0.0:
            share_deposits_invested = np.nan
        else:
            share_deposits_invested = \
                total_investments / self.assets["deposits"][self.firm_bank_id]
            assert 0.0 <= share_deposits_invested <= 1.0, pdb.set_trace()

        self.share_deposits_invested.append(share_deposits_invested)
        self.assets["deposits"][self.firm_bank_id] -= total_investments

        assert self.assets["deposits"][self.firm_bank_id] >= 0, \
            "Deposits after investments should be >= 0, not {}".format(
                self.assets["deposits"])

        innovation_success, spillover_success, information_success = \
            self.compute_success(
                rd_investment=rd_investment,
                info_investment=info_investment,
                spillover_investment=spillover_investment,
                target_node=target_node)

        # 9. Go to new spot on the product space
        new_product, updated_innovation_capabilities, \
        updated_spillover_capabilities = \
            self.make_move(
                target_node=target_node,
                innovation_success=innovation_success,
                spillover_success=spillover_success,
                complexity=target_complexity,
                firms=target_firms,
                same_node=same_node)

        # 10. Update capital stock and account
        updated_capital_stock, updated_account_firm = self.update_capital_stock(
            k_investment=k_investment,
            payback=payback,
            credit=credit_granted,
            profit=self.profit[-1],
            remaining_profit=remaining_profit)

        # 11. Update firm properties
        new_range_of_vision = self.update_range_of_vision(
            information_success=information_success)
        if self.print_progress:
            print(self, self.model_instance.t, ": ", new_range_of_vision)

        self.k_investment = k_investment
        self.rd_investment = rd_investment
        self.ac_investment = spillover_investment
        self.innovation_capabilities.append(updated_innovation_capabilities)
        self.spillover_capabilities.append(updated_spillover_capabilities)
        self.capital_stock.append(updated_capital_stock)
        self.account_firm.append(updated_account_firm)
        self.range_of_vision = new_range_of_vision
        if self.product_space_position == new_product:
            self.remained_same_node = True
        else:
            self.remained_same_node = False
        self.product_space_position = new_product
        self.path.append(new_product)
        firm_balance_sheet = self.get_balance_sheet()
        self.firm_deposits = firm_balance_sheet["deposits"]
        self.firm_net_worth = firm_balance_sheet["net_worth"]
        self.firm_assets = firm_balance_sheet["assets"]
        self.firm_liabilities = firm_balance_sheet["liabilities"]

        if self.ultimate_target_node == self.product_space_position:
            self.ultimate_target_node = False

        if self.print_progress == True:
            print(file_marker + self.path)

        self.check_final_spot()

        return updated_capital_stock, updated_account_firm

    def choose_target_node(self, price, output):
        """Find the target node for the firm

        Firm chooses from the nodes it sees the one that promises highest
        expected profits as ultimate target. This is the 
        'ultimate target node'.
        However, since a firm can only go one step at a time, it cannot
        reach ultimate target nodes that are more than one neighbor 
        away. In this case, the firm needs to 'save moves': it remains
        on its current product and in the following round it has 
        the possibility to move one step further. Once it has saved 
        enough moves, it goes directly to the ultimate target node. 
        This way, we prevent firms from being stuck on unattractive 
        intermediate products (or changing their production more often 
        than necessary).
        
        The overall procedure is as follows: 
            1. Identify the ultimate target node
            2. Check the shortest path to the ultimate target node
            3. Get the immediate next step (i.e. product to move next)
            4. Collect information about the target node and return it
        
        Parameters
        ----------
        price : float
            price of the current product
        output : float
            current output of the firm

        Returns
        -------
        target_node : int
            the next node that a firm aims to go to
        same_node : bool
            Does the firm wish to stay in its current position?
        target_complexity : float
            The complexity of the target node
        target_firms : int
            The number of firms on the target node
        """
        # 1. Identify the ultimate target node
        if self.ultimate_target_node is False:
            # The firm has not yet selected a target node
            """Gather information about price, complexity and output
            of all products in sight"""
            price_dict = {}
            output_dict = {}
            complexity_dict = {}

            price_dict[self.product_space_position] = price
            output_dict[self.product_space_position] = output
            complexity_dict[self.product_space_position] = \
                self.model_instance.product_space.nodes[ \
                    self.product_space_position]["complexity"]

            visible_products = self.model_instance.get_all_visible_products(
                self.product_space_position, self.range_of_vision)

            if self.cant_go_back == True:
                for p in range(len(visible_products)):
                    if visible_products[p] in self.path:
                        visible_products.pop(p)

            for prod in visible_products:
                price_dict[prod] = \
                    self.model_instance.product_space.nodes[prod]["price"]
                output_dict[prod] = min(
                    self.model_instance.product_space.nodes[prod]["capital"]
                    * self.capital_productivity,
                    self.nominal_demand)
                complexity_dict[prod] = \
                    self.model_instance.product_space.nodes[prod]["complexity"]

            """Compute the expected profits based on the info"""
            expected_profit_dict = self.compute_expected_profits(
                products_of_interest=visible_products,
                price_dict=price_dict,
                output_dict=output_dict)

            """Pick product with maximum expected profit"""
            max_expected_profit = max(expected_profit_dict,
                                      key=expected_profit_dict.get)
            ultimate_target_node = max_expected_profit

        else:
            # The firm has selected a target node that so far it was
            # unable to reach.
            ultimate_target_node = self.ultimate_target_node

        # 2. Specify the path to the ultimate target node
        path = nx.shortest_path(self.model_instance.product_space,
                                source=self.product_space_position,
                                target=ultimate_target_node)

        # 3. Get the immediate next step (i.e. product to move next)
        if len(path) == 1:
            # The firm wishes to stay on current position
            target_node = path[0]
            same_node = True
            self.wait_for_move = False
            self.available_steps = 1
            assert target_node == self.product_space_position, \
                "Internal error in choose_target_node"
        elif len(path) == 2:
            # The firm wishes to move to a neighbor product
            target_node = path[1]
            same_node = False
            self.wait_for_move = False
            self.available_steps = 1
        else:
            ultimate_target_node_distance = len(path) - 1
            if self.available_steps >= ultimate_target_node_distance:
                target_node = ultimate_target_node
                if target_node == self.product_space_position:
                    same_node = True
                else:
                    same_node = False
                self.wait_for_move = False
            else:
                self.ultimate_target_node = ultimate_target_node
                target_node = path[0]
                same_node = True
                self.wait_for_move = True
                self.available_steps += 1

        # 4. Collect information about the target node and return it
        target_complexity = \
            self.model_instance.product_space.nodes[target_node].get(
                "complexity")
        target_firms = \
            self.model_instance.product_space.nodes[target_node].get(
                "firms")

        return target_node, same_node, target_complexity, target_firms

    def compute_expected_profits(self, products_of_interest, price_dict,
                                 output_dict):
        """Firms compute profits for all visible products.
        
        1. Gather price and total output
        2. Compute expected profits given the current price and output
        3. Get the complexity of the target & its current market size
        4. Compute the desired investment if the firm were on the new spot

        Parameters
        ----------
        products_of_interest : dict
            A dictionary with products as keys and their distances as 
            values. Usually the set of products that a firm can see at 
            the moment of choosing a new location on the product space.
        price_dict : dict
            dict with info on the prices of the `products_of_interest`
        output_dict : dict
            dict with info on production of the `products_of_interest`

        Returns
        -------
        expected_profit_dict : dict
            dictionary that gives information on the expected profit of 
            each visible product
        """
        expected_profit_dict = {}
        for product in products_of_interest:
            # 1. Gather price and total output 
            target_price = price_dict.get(product)
            target_total_output = output_dict.get(product)

            if isinstance(target_price, (np.float64, float)) == False:
                pdb.set_trace()

            if target_total_output == 0:
                market_share = 1
            else:
                market_share = \
                    self.capital_stock[-1] / (
                            self.model_instance.product_space.nodes[
                                product]["capital"] + self.capital_stock[-1])
            output = max(min(self.nominal_demand * market_share,
                             self.capital_stock[
                                 -1] * self.capital_productivity),
                         0)

            # 2. Compute expected profits given the current price and output
            expected_price_cost_ratio, expected_profit, \
            expected_payback, expected_gross_revenue, expected_interest = \
                self.get_profit(
                    price=target_price,
                    output=output)

            expected_profit_dict[product] = expected_profit

        return expected_profit_dict

    def update_range_of_vision(self, information_success):
        """Computes new range of vision

        Depending on whether absorptive capacities were successful, 
        the range of vision is broadened such that the firm can see one 
        further product.

        Parameters
        ----------
        information_success : boolean
            gives information on whether absorptive capacities that make 
            it easier to evaluate the environment where successful

        Returns
        -------
        range_of_vision : int
            gives info on how many products a firm has information on
        """
        if information_success == True:
            range_of_vision = min(
                self.max_info * len(self.model_instance.product_space),
                self.range_of_vision + 0.5)
        else:
            range_of_vision = self.range_of_vision

        return range_of_vision

    def make_move(self, target_node, innovation_success, spillover_success,
                  complexity, firms, same_node):
        """Can the target node be reached?

        Successful capability measures will lead to higher capabilities 
        in the respective categories.
        The firm can make a move away from its current position if
        (1) R&D capabilities are sufficient or if
        (2) spillover capabilities are sufficient AND there is already a
        market for the target product (i.e. there are firms that can
        spillover their information)

        Parameters
        ----------
        target_node : int
            the node that the firm would like to go to
        innovation_success : bool
            gives info on whether R&D investment was successful
        spillover_success : bool
            gives info on whether absorptive capacities, which make 
            spillovers easier, were successful
        firms : int
            the number of firms that produce the target product
        complexity : float
            the complexity of the target product
        same_node : bool
            does the firm wish to stay in its current market?

        Returns
        -------
        new_product : int
            the product that the firm will produce in the next period

        innovation_capabilities : float
            the new capability level associated with R&D

        spillover_capabilities : float
            the new capability level associated with absorptive capacities
        """
        if innovation_success:
            innovation_capabilities = self.innovation_capabilities[-1] + \
                                      self.effect_rd_success

        else:
            innovation_capabilities = self.innovation_capabilities[-1]

        if spillover_success:
            spillover_capabilities = self.spillover_capabilities[-1] + \
                                     self.effect_spillover_success
        else:
            spillover_capabilities = self.spillover_capabilities[-1]

        if self.wait_for_move == False and same_node == False:
            if complexity <= innovation_capabilities:
                new_product = target_node
                self.available_steps = 1
            elif firms > 0 and complexity <= spillover_capabilities:
                new_product = target_node
                self.available_steps = 1
            else:
                new_product = self.product_space_position
        else:
            new_product = self.product_space_position

        return new_product, innovation_capabilities, spillover_capabilities

    def output_and_capability_costs(self, previous_capital_stock, price):
        """Determines current output
        
        First, the market share of the firm gets computed. 
        Usually, the market share is computed as the total capital stock
        of the firm divided by the aggregated capital stock of all firms
        that are located on this position on the product space (in order
        to account for firm size).
        If the product has not been produced in the previous period at 
        all, the market share is set to one.
        Second, the actual output gets determined. To this end, the firm
        first checks whether the production is worth it by testing 
        whether the price is at least equal to the quotient of the cost
        share and capital productivity.
        If it is profitable, it produces its (capital stock * capital
        productivity - inventories). If total production would
        exceed the maximum total demand, then the firm produces
        according to its market share, taking inventories into account.

        Parameters
        ----------
        previous_capital_stock : float
            The firm's capital stock at the end of the last period.
        price : float
            The current price of the product the firm is producing

        Returns
        -------
        output : float
            The firm's output in the current time step.
        innovative_rd : float
        """

        """1. Compute market share """
        if len(self.total_output_prod_produced) == 0:
            market_share = \
                1 / self.model_instance.product_space.nodes[
                    self.product_space_position].get("firms")
        else:
            if self.total_output_prod_produced[-1] == 0:
                market_share = 1
            else:
                market_share = self.capital_stock[-1] / \
                               (self.model_instance.product_space.nodes[
                                    self.product_space_position]["capital"] +
                                self.capital_stock[-1])
        """2. Determine output"""
        if price < self.cost_share / self.capital_productivity:
            output = 0
        else:
            output = max(
                min(self.capital_productivity * previous_capital_stock -
                    self.inventories[-1],
                    self.nominal_demand * market_share - self.inventories[-1]
                    ),
                0)

        if self.print_progress == True:
            print(file_marker + 'Firm ', self.id, 'output ', output)

        return output

    def get_profit(self, price, output):
        """Compute profits for given prices and output values
        
        1. Compute the amount of interest the firm has to pay to the 
            bank, and the the amount of money it has to pay as payback 
            to the bank
        2. Compute profit as the difference between gross revenue and
            interest and installment payments
        3. Compute the price-cost-ratio

        Parameters
        ----------
        price : float
            Current market price.

        output : float
            The firm's output in the current step.

        Returns
        -------
        price_cost_ratio : float
            Will later be used to determine a firm's investment.
        profit : float
            The firm's profit.
        payback : float
            The sum of loans that have to be paid back to the bank.
        gross_revenue : float 
            Price times output, i.e. the money the firm makes by selling
            its output
        interest_payment : float
            The amount of interest on loans.
        """
        if self.account_firm[-1] >= 0:
            interest_payment = self.account_firm[-1] * self.rate_on_deposits * (
                -1)
            payback = 0
            assert interest_payment <= 0, \
                "Firm account positive, but interest payment enters negatively"
        else:
            interest_payment = -(self.account_firm[-1] * self.rate_on_loans)
            # this will be a positive number
            assert interest_payment >= 0, \
                "Firm account negative, but interest payment enters positively"
            payback = -self.account_firm[-1] * self.payback_rate
            assert payback > 0, "Payback is negative number: {}".format(payback)

        assert (price * output) >= 0, "Price times output is negative!"
        assert (self.cost_share * output / self.capital_productivity) >= 0, \
            "Cost share is a negative number!"

        gross_revenue = (price * output) - \
                        (self.cost_share * output / self.capital_productivity)

        profit = gross_revenue - interest_payment - payback

        price_cost_ratio = price * self.capital_productivity / self.cost_share

        return price_cost_ratio, profit, payback, gross_revenue, interest_payment

    def get_demand_for_investment(self, output, deposits, target_complexity,
                                  target_firms):
        """Computes the desired investment of the firm
        
        1. Compute demand for capital investment
        First, compute `demand_option` as the potential to increase 
        production for the product the firm is producing.
        Second, compute the `capital gap`. To this end, take the previous 
        output plus the `demand_option` minus the inventories - i.e. the 
        total  to be produced output - divided by the capital 
        productivity. This gives the amount of capital necessary to
        produce what is to be produced. From this one subtracts the
        amount of capital still present after depreciation. 
        The demanded investment into capital is the capital gap, or,
        in case the firm has selected an ultimate target node,
        the maximum of the current capital gap and the depreciated capital
        stock (this is to reflect the uncertainty that a firm faces before
        changing markets).
        
        2. Compute demand for R&D und spillover investment
        - When the firm has not yet selected an ultimate target node, d
            emand for R&D and spillover will be 0. This is to reflect 
            that firms that do not wish to change markets have no need 
            to invest in new capabilities.
        - When the target complexity is smaller or equal the current
            innovation or spillover capabilities, there is also no demand
            for R&D or spillovers. This is because, in this case, the firm
            already knows how to produce the target product and, therefore,
            does not need to invest into any further capabilities for now.
        - In all other cases, the demand for R&D is the maximum of the 
            `target_prob´('self.p_information'), deposits
            times self.p_innovation divided by the sum of 
            self.p_innovation + self.p_spillovers + self.p_information. 
            For the demanded spillovers its the same, only with 
            self.p_spillovers instead of self.p_innovation.
            Finally, two adjustments are potentially made:
            if there are no target firms, demanded spillovers are zero.
            If there are more than 10 target firms, demanded R&D is zero.
            Then R&D and Spillovers are adjusted according to the nb of 
            firms.
        
        Parameters
        ----------
        output : float
            The output of the firm in the current time step
        deposits : float
            The deposits of the firm in the current period
        target_complexity : float
            The complexity of the target node
        target_firms : int
            The number of firms on the target node
        
        Returns
        -------
        demand_investment : float
            Total investment demand
        demand_credit : float
            The firm's demand for credit
        demand_information : float
            The firm's demand for investment into its range of vision
        demand_spillovers : float
            The firm's demand for investment into spillovers
        demand_k : float
            The firm's demand for investment into its capital stock
        demand_rd : float
            The firm's demand for investment into R&D
        """

        """Compute demand for capital investment"""
        demand_option = max(
            0,
            (self.nominal_demand - self.total_output_prod_produced[-1])
        )
        capital_gap = (output + demand_option - self.inventories[-1]) / \
                      self.capital_productivity - \
                      (1 - self.depreciation_rate) * self.capital_stock[-1]
        if self.ultimate_target_node == False:
            demand_k = max(0, capital_gap)
        else:
            demand_k = max(
                self.depreciation_rate * self.capital_stock[-1],
                capital_gap, 0)

        """Determine amount of R&D and spillovers demanded"""
        if target_complexity == None and target_firms == None:
            demand_rd = 0
            demand_spillovers = 0
        elif target_complexity <= self.innovation_capabilities[-1] or \
                target_complexity <= self.spillover_capabilities[-1]:
            demand_rd = 0
            demand_spillovers = 0
        else:
            demand_rd = max(
                self.target_prob / self.p_innovation,
                deposits * self.p_innovation /
                (self.p_innovation + self.p_spillovers + self.p_information)
            )
            demand_spillovers = max(
                self.target_prob / self.p_spillovers,
                deposits * self.p_spillovers /
                (self.p_innovation + self.p_spillovers + self.p_information)
            )
            if target_firms == 0:
                demand_spillovers = 0
            elif target_firms > 10:
                demand_rd = 0
            else:
                demand_rd = demand_rd * (1 - (target_firms / (
                        10 + 1)))
                demand_spillovers = demand_spillovers * target_firms / 10

        demand_information = max(
            self.target_prob / self.p_information,
            profit * self.p_information /
            (self.p_innovation + self.p_spillovers + self.p_information))

        """Determine amount of investment demanded"""
        demand_investment = \
            demand_k + demand_rd + demand_information + demand_spillovers

        if self.print_progress and demand_investment <= profit:
            print("only want to invest {} per cent of deposits.".format(
                demand_investment / profit))

        """Determine amount of credit demanded"""
        demand_credit = max(0.0, demand_investment - profit)
        self.demand_rd.append(demand_rd)

        return demand_investment, demand_credit, demand_information, \
            demand_spillovers, demand_k, demand_rd

    def get_investment(self, financial_constraint, demand_k, demand_rd,
                       demand_info, demand_spillovers):
        """Computes actual investment

        Here, the demand for investment is aligned with the firm's financial
        means.
        (1) The investment into its capital stock has highest priority,
        i.e. the firm will invest as much as is financially possible 
        (and as is demanded) into capital (indicating a higher priority 
        of production than capability accumulation)
        (2) Then the financially possible investment into capability 
        measures is computed. If total capability-investment remain under 
        the financial constraint, the entire demand will be realized. 
        Otherwise, actual investment will be split up between the 
        measures, relative to the size of their demand.

        Parameters
        ----------
        financial_constraint : float
            Amount of money the firm can spend in this period
        demand_k : float
            The firm's demand for investment into its capital stock
        demand_rd : float
            The firm's demand for investment into R&D
        demand_info : float
            The firm's demand for investment into its range of vision
        demand_spillovers : float
            The firm's demand for investment into spillovers

        Returns
        -------
        rd_investment : float
            The firm's actual investment into R&D
        info_investment : float
            The firm's actual investment in its range of vision
        spillover_investment : float
            The firm's actual investment in spillovers
        k_investment : float
            The firm's actual investment in its capital stock
        remaining_profit : float
            Profit after investment
        """
        k_investment = min(demand_k, financial_constraint)

        total_capability_demand = \
            demand_info + demand_rd + demand_spillovers
        capability_constraint = financial_constraint - k_investment
        if self.print_progress:
            print("financial_constraint - total investment: {}".format(
                financial_constraint - (
                        total_capability_demand + k_investment))
            )
        info_investment = min(
            demand_info,
            (demand_info / total_capability_demand) * capability_constraint)
        rd_investment = min(
            demand_rd,
            (demand_rd / total_capability_demand) * capability_constraint)
        spillover_investment = min(
            demand_spillovers,
            (1 - (demand_rd + demand_info) / total_capability_demand) * \
            capability_constraint)

        total_investments = \
            k_investment + info_investment + rd_investment + spillover_investment

        assert total_investments - financial_constraint < 0.00000005, \
            "firms invest more than they can"
        remaining_profit = max(0.0, financial_constraint - total_investments)

        return rd_investment, info_investment, spillover_investment, \
            k_investment, remaining_profit

    def update_capital_stock(self, k_investment, payback, credit, profit,
                             remaining_profit):
        """Updates the firm's capital stock and bank account
                
        Parameters
        ----------
        k_investment : float
            Investment into the capital stock
        payback : float
            The money the firm has paid to the bank in order to pay
            back debts. It is already included in `profit`: when 
            computing the latter in `get_profit()`, the `payback` has 
            been subtracted from gross revenues.
        credit : float or int
            Amount of credit that was granted to the firm
        profit : float
            The firm's profit
        remaining_profit : float
            The firm's profit after investment
        
        Returns
        -------
        updated_capital_stock : float
            The firm's capital stock after this period's depreciation and
            investment into new capital
        updated_account_firm : float
            The firm's overall account after all transactions in this period
            have been executed
        """
        updated_capital_stock = \
            (1 - self.depreciation_rate) * self.capital_stock[-1] + k_investment
        assert updated_capital_stock > 0, \
            "capital stock became negative!"

        if k_investment == 0:
            assert updated_capital_stock < self.capital_stock[-1]

        self.assets["capital_stock"] = updated_capital_stock

        added_loss = 0
        if profit < 0:
            added_loss = profit

        # payback: Positive nb, increases account
        # added_loss: Negative number, decreases account
        # credit: positive number, increases account???
        # remaining_profit: positive number, increases account
        assert payback >= 0, \
            "payback is a negative nb: {}".format(payback)
        assert added_loss <= 0, \
            "added_loss is a positive nb: {}".format(added_loss)
        assert credit >= 0, \
            "Credit is a negative nb: {}".format(credit)
        assert remaining_profit >= 0, \
            "remaining_profit is negative: {}".format(remaining_profit)
        updated_account_firm = \
            self.account_firm[-1] + payback + added_loss - credit + \
            remaining_profit

        return updated_capital_stock, updated_account_firm

    def compute_success(self, rd_investment, info_investment,
                        spillover_investment, target_node):
        """Conduct capability accumulation activity
        
        Firms now
        (1) Conduct R&D. The higher their investment, the higher their
        chances of a success.
        (2) Take measures to extend their absorptive capacities:
            a) They broaden their information on the product space (range of
            vision)
            b) They learn from other firms (spillovers)

        Parameters
        ----------
        rd_investment : float
            The firm's actual investment into R&D
        info_investment : float
            The firm's actual investment in its range of vision
        spillover_investment : float
            The firm's actual investment in spillovers
        target_node : int
            The firms ultimate target node

        Returns
        -------
        innovation_success : bool
            True if R&D was successful
        spillover_success : bool
            True if spillover-measures were successful
        information_success : bool
            True if extension of range of vision was successful
        """
        # 1st step: was R&D successful?
        if rd_investment < 1 / self.p_innovation:
            success_prob_innovation = 0
        else:
            success_prob_innovation = self.p_innovation - 1 / rd_investment

        assert (success_prob_innovation <= 1 and success_prob_innovation >= 0), \
            "Computation of success probability yields figure l" \
            "arger than 1 or less than 0"

        if success_prob_innovation > 1:
            success_prob_innovation = 1
        elif success_prob_innovation < 0:
            success_prob_innovation = 0

        if bernoulli.rvs(success_prob_innovation, size=1) == 1:
            innovation_success = True
        else:
            innovation_success = False

        # 2nd step: was trying to acquire spillovers successfull?
        if spillover_investment < 1 / self.p_spillovers:
            success_prob_spillovers = 0
        else:
            success_prob_spillovers = \
                self.p_spillovers - 1 / spillover_investment
        assert success_prob_spillovers <= 1,\
            "Success probability yields figure larger than 1."

        assert success_prob_spillovers >= 0, \
            "Success probability yields figure less than 0."

        neighboring_firms = self.model_instance.product_space.nodes[
            target_node].get("firms")

        spillover_success = False
        for n in range(neighboring_firms):
            if bernoulli.rvs(success_prob_spillovers, size=1) == 1:
                spillover_success = True
                break
            else:
                continue

        # 3rd step: was trying to broaden information on the network
        # successful?
        information_success = False
        if info_investment < 1 / self.p_information:
            success_prob_information = 0
        else:
            success_prob_information = self.p_information - 1 / info_investment
        if success_prob_information > 1:
            success_prob_information = 1
        elif success_prob_information < 0:
            success_prob_information = 0
        if bernoulli.rvs(success_prob_information, size=1) == 1:
            information_success = True

        return innovation_success, spillover_success, information_success
