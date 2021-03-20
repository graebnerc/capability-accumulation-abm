import os
file_marker = "[" + str(os.path.basename(__file__)) + "]: "

class Bank:
    """Banks lend money to firms.
    
    Firms put money they have left to the bank (deposits), but also ask 
    the bank for credit if they want to invest more than they earn.
    
    The bank pays interest on deposits and earns interest for loads.
    
    Updating the bank means to update total amount of loans and deposits.
    """
    def __init__(self, financial_regime_parameter, rate_on_deposits,
                 rate_on_loans, payback_rate, identifier, model_instance):
        self.model_instance = model_instance
        self.identifier = identifier
        self.financial_regime_parameter = financial_regime_parameter
        self.account_bank = [0]
        self.rate_on_deposits = rate_on_deposits
        self.rate_on_loans = rate_on_loans
        self.payback_rate = payback_rate
        self.assets = {
            "loans": {},
            "deposits": 0.0
        }
        self.liabilities = {
            "deposits": {}
        }
        self.bank_net_worth = 0.0
        self.bank_assets = 0.0
        self.bank_liabilities = 0.0
        
    def get_balance_sheet(self):
        """Returns the balance sheet of the bank
        
        Computes net worth as assets minus liabilities.
        
        Returns
        -------
        dict
            Dictionary with information on `net_worth`, `assets`,
            and `liabilities`.
        """
        total_assets = self.assets["deposits"] + \
            sum(list(self.assets["loans"].values()))
        assert total_assets >= 0, "Total assets of bank are negative"
        total_liabilities = sum(self.liabilities["deposits"].values())
        assert total_liabilities <= 0, \
            "Total liabilities of bank are negative"
        net_worth = total_assets + total_liabilities
        balance_sheet = {
            "net_worth": net_worth,
            "assets": total_assets,
            "liabilities": total_liabilities
        }
        return balance_sheet
        
    def get_interest_on_deposits(self, firm_id):
        """Computes interest a firm gets on its deposits

        Parameters
        ----------
        firm_id : int
            id of the firm that holds deposits
        
        Returns
        -------
        float
            The total amount of interests the firms receives
        """
        amount_deposits = abs(self.liabilities["deposits"][firm_id])
        interest_payment = amount_deposits * self.rate_on_deposits
        return interest_payment  
     
    def get_interest_for_loans(self, firm_id):
        """Gets amount firm has to pay interest on its loans

        Parameters
        ----------
        firm_id : int
            id of the firm that holds loans
        
        Returns
        -------
        float
            The total amount of interests the firms has to pay
        """
        amount_loan = abs(self.assets["loans"][firm_id])
        interest_payment = amount_loan * self.rate_on_loans

        assert interest_payment >= 0, "Firms have to pay for loans"
        return interest_payment       
        
    def receive_loan_payback(self, firm_id, payback_amount):
        """Receives payback for existing loan

        Parameters
        ----------
        firm_id : int
            The id of the firm that pays back the loan
        payback_amount : float
            The amount the firm pays back
        
        Raises
        ------
        AssertionError
            If loans of bank have wrong sign
        """
        if payback_amount > self.assets["loans"][firm_id]:
            payback_amount = self.assets["loans"][firm_id]
        self.assets["deposits"] += payback_amount
        self.assets["loans"][firm_id] -= payback_amount
        new_loan_level = self.assets["loans"][firm_id]

        assert new_loan_level >= 0, \
            "Loans of bank must not be smaller than zero: {}".format(
                new_loan_level)

    def update_bank(self, all_deposits, all_loans):
        """Update all relevant parameters for the single banks

        Bank account is updated here at the end of each time step

        Parameters
        ----------
        all_deposits : float
            Sum of all the deposits that the firms have with the bank
        all_loans : float
            Sum of all loans given to the banks

        Returns
        -------
        updated_account_bank : float
            The bank's account after all transactions have been executed

        """        
        interest_payment = all_deposits * self.rate_on_deposits + \
            all_loans*self.rate_on_loans

        updated_account_bank = all_deposits + all_loans

        self.account_bank.append(updated_account_bank)

        return updated_account_bank

    def credibility(self, firm_id, profit, capital_stock,
                    credit_demand, deposits):
        """Determines credit that is granted to the firm
        
        The maximal amount of credit (`maximum_credit`) is given by the
        product of the firm's profit rate and the 
        `financial_regime_parameter`.
        It the demanded credit exceeds the maximum credit, a credit of 
        the amount `maximum_credit` is granted. Otherwise the demanded 
        sum is granted. If, for some case, the demanded sum is negative 
        then an error is raised.
        
        In case the credit is granted, it gets already written into 
        the balance sheet of the firm and bank.
        
        Parameters
        ----------
        firm_id : int
            The id of the firm that applies for a loan
        profit : float
            The profits of the firm that applies for credit.
        capital_stock : float
            The capital stock of the firm
        credit_demand : float
            The amount of credit demanded by the firm.
        
        Returns
        -------
        float
            The credit provided by the bank
        
        Raises
        ------
        AssertionError
            If `credit_demand` is smaller than zero.
        """
        assert credit_demand >= 0.0, \
            "Demanded credit is negative: {}.".format(credit_demand)
        rate_of_return = profit / capital_stock
        maximum_credit = self.financial_regime_parameter * deposits

        if credit_demand >= maximum_credit:
            credit = max(0, maximum_credit)
        elif profit < 0:
            credit = 0
        else:
            credit = credit_demand
        
        if credit > 0:
            self.model_instance.place_loan(firm_id, self.identifier, credit)

        return credit
