'''
cbc_hb.py
-----------
This library of functions does choice based conjoint analysis using 
heirarchical bayes as implemented by Sawtooth Software.
Author: William Ma 
Date: Feburary 28, 2018
References:
    https://www.sawtoothsoftware.com/download/techpap/CBCHB_Manual.pdf
'''
import logging
import sys

import pandas as pd
import numpy as np

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

SPECIAL_STATE = np.random.RandomState()
np.seterr(all='raise')

class Cbc_hb:
    """Create model to run choice based heirarchical bayes conjoint analysis
    
    Class Variables
    ---------------
    seed_state : tuple
        Tuple containing seed state for numpy
    
    Instance Variables
    ------------------
    alphas : Numpy array
        Numpy array containing the mean partworth of parameters and (the first
        one) the normalizing term
        Dimension are (iters x n_parameters-n_attributes+1)
    attribute_levels : list
        List of list of the parameter levels grouped by attribute
    attribute_names : list
        List of the attribute names being tested
    avg_betas : Numpy array
        Average of betas drawn
        Dimensions are (iters x n_parameters-n_attributes+1)
    betas : Numpy array
        Numpy array containing the partworths of parameters for each respondent
        and (the first one) the normalizing term
        Dimensions are (iters x n_parameters-n_attributes+1 x n_respondents)
    choices : Numpy array
        Matrix of the choices that each respondent made for each comparison
        Dimensions are (n_comparisons x n_respondents x n_parameters)
    comparisons : Pandas Dataframe
        Dataframe containing the comparisons data
    iters : int
        Number of iterations to run the model
    n_attributes : int
        Number of attributes
    n_comparisons: int
        Number of comparisons
    n_options : int
        Number of choices in each comparison
    n_parameters : int
        Total number of parameters
    n_attr_parameters : Numpy array
        Vector of the number of parameters per attribute
        Dimensions are (n_attributes)
    n_respondents : int
        Number of respondents
    options : Numpy array
        Matrix containing the binary of each choice for every comparison
        Dimensions are (n_options x n_comparisons x n_parameters)
    parameter_names : list
        List of the parameter names
    probabilities : Numpy array
        Vector of the probabilities based off of the avg_betas
        Dimensions are (n_parameters)
    profiles : Pandas DataFrame
        Dataframe containing the profiles presented
    profile_bin : Pandas DataFrame
        Dataframe containing the binary data of the profiles
    relative_importance : Numpy array
        Vector containing the relative importance of each attribute
        Dimensions are (n_attribues)
    selections : Pandas DataFrame
        Dataframe containing the selections made for each comparison for each
        respondent
    _choices : Numpy matrix
        Modified version of choices that removes one parameter of each 
        attribute for normalization and adds one column of ones for the 
        normalizing term
        Dimensions are (n_comparisons x n_respondents x n_parameters-n_attributes+1)
    _removed : list
        List of the indicies that were removed from _choices and _options
        Dimensions are (n_attributes)
    _options : Numpy matrix
        Modified version of options that removes one parameter of each 
        attribute for normalization and adds one column of ones for the 
        normalizing term
        Dimensions are (n_options x n_comparisons x n_parameters-n_attributes+1)

    Parameters
    ----------

    Returns
    -------

    """

    seed_state = SPECIAL_STATE.get_state()

    def __init__(self, profiles, comparisons, selections):
        '''
        Create object that contains all the data to run conjoint analysis

        Parameters
        ----------
        profiles : DataFrame
            Dataframe with column names of the variations, first column is 
            profile number, each of the columns contains the variation used
        comparisons : DataFrame
            DataFrame with the following columns with headers:
                * First column is the comparison number
                * Remaining columns are the choices that were presented in comparison
        selections : DataFrame
            DataFrame with the following columns:
                * First column is the comparison number
                * Next n_respondents columns is the choice of the i-th individual
        '''
        # Save arguments to corresponding instance variables
        self.profiles = profiles
        self.comparisons = comparisons
        self.selections = selections 

        # Convert the profiles into binary values
        self.profile_bin = pd.get_dummies(profiles)#, drop_first=True)
        
        # Save parameter names from profile_bin
        self.parameter_names = self.profile_bin.columns

        # Get attribute names and parameter counts
        self.attribute_names = list()
        self.attribute_levels = list()
        self.n_attr_parameters = np.array([])
        for param_str in self.parameter_names:
            tokens = param_str.split("_")
            param_name = tokens[0]
            param_level = tokens[1]
            if param_name not in self.attribute_names:
                self.attribute_names.append(param_name)
                self.attribute_levels.append(list())
                self.n_attr_parameters = np.append(self.n_attr_parameters, 0)
            ind = self.attribute_names.index(param_name)
            self.n_attr_parameters[ind] += 1
            self.attribute_levels[ind].append(param_level)

        self.n_parameters = int(sum(self.n_attr_parameters))
        self.n_attributes = len(self.attribute_names)
        self.n_respondents = len(selections.columns)  # Respondants
        self.n_comparisons = comparisons.shape[0]     # Number of comparisons
        self.n_options = len(comparisons.columns)     # Number of options in 
                                                      # each comparison

        # Convert the choices into binary values
        profile_options = dict()
        for profile in comparisons.columns:
            profile_values = self.profile_bin.loc[comparisons[profile]]
            profile_options[profile] = profile_values.reset_index(drop=True)
        self.options = pd.concat(profile_options, axis=1)
        self.options = np.stack([self.options['Profile{}'.format(i+1)] 
                                    for i in range(self.n_options)]) 

        # Create choices 3D array
        choice_arrs = list()
        for comparison in range(comparisons.shape[0]):
            comp_arrs = list()
            for selection in selections.iloc[comparison, :].values:
                if selection >= 0:
                    comp_arrs.append(self.options[selection, comparison, :])
                else:
                    comp_arrs.append(np.zeros(self.n_parameters))
            choice_arrs.append(np.stack(comp_arrs))
        self.choices = np.stack(choice_arrs)

        # TODO: Use np.vstack instead of appending and deleting arrays
        self._choices = self.choices
        self._options = self.options
        self._removed = list()
        ind = 0
        for i in range(self.n_attributes):
            # Subtract by i to account for the columns already removed
            self._options = np.delete(self._options, (ind-i), axis=2)
            self._choices = np.delete(self._choices, (ind-i), axis=2)
            self._removed.append(ind)
            ind += self.n_attr_parameters[i]
        options_first = np.ones(self.n_comparisons)
        choices_first = (np.ones(self.n_comparisons * self.n_respondents)
                            .reshape(self.n_comparisons, self.n_respondents))
        self._options = np.insert(self._options, 0, options_first, axis=2)
        self._choices = np.insert(self._choices, 0, choices_first, axis=2)

        # Initialize remaining instance variables to zero or empty
        self.iters = 0
        self.alphas = np.array([])
        self.betas = np.array([])
        self.avg_betas = np.array([])
        self.probabilities = np.zeros(self.n_parameters)
        self.relative_importance = np.zeros(self.n_attributes)


    def model(self, iters=20000, proportionality_factor=0.1, acceptance_rate=0.3):
        """Takes in data and does conjoint analysis on it
            
        Parameters
        ----------
        iters : int
            The number of iterations to run the simulation for (Default value = 20000)
        proportionality_factor :
            Proportionality or jump factor to mupltiply by (Default value = 0.1)
        acceptance_rate :
            Proportion of betas we want to accept (Default value = 0.3)

        Yields
        ------
        self.alphas : Numpy array
            Numpy array containing the mean partworth of included parameters
            and one for the normalizing term
            Dimensions are (iters, n_parameters-n_attributes+1)
        self.betas : Numpy array
            Numpy array containing the partworths of included parameters for
            each respondent and one for the normalizing term
            Dimensions are (iters, n_parameters-n_attributes+1, n_respondents)
        self.avg_betas : Numpy array
            Numpy array containing the mean partworth of included parameters
            and one for the normalizing term
            Dimensions are (iters, n_parameters-n_attributes+1)

        Returns
        -------

        Notes
        -----
        [1] https://www.sawtoothsoftware.com/download/techpap/CBCHB_Manual.pdf
        """
        
        # Define instance variable of iters
        self.iters = iters

        # Initialize Numpy arrays to yield
        self.alphas = np.zeros((iters, self.n_parameters-self.n_attributes+1))
        self.betas = np.zeros((iters, self.n_parameters-self.n_attributes+1, self.n_respondents))
        self.avg_betas = np.zeros((iters, self.n_parameters-self.n_attributes+1))

        # We use an all the priors are zeros
        alpha = np.zeros(self.n_parameters-self.n_attributes+1)
        beta = np.zeros((self.n_parameters-self.n_attributes+1, self.n_respondents))
        D = np.eye(self.n_parameters-self.n_attributes+1)

        for iter_num in range(iters):
            logging.debug("\n********************* ITER {} ****************\n".format(iter_num))

            avg_beta = np.mean(beta, axis=1)
            cov = (1.0/self.n_respondents) * D
            alpha = SPECIAL_STATE.multivariate_normal(np.transpose(avg_beta), cov)

            # logging.debug("alpha_avg: {}".format(avg_beta))
            # logging.debug("alpha_cov: {}".format(cov))
            # logging.debug("alpha: {}".format(alpha))

            # N = respondents + (n_parameters - n_attributes)
            big_N = self.n_respondents + \
                    self.n_parameters - self.n_attributes + 1

            # H = pI + sum_{i=0}^{n_respondents} (alpha - beta_i)(alpha - beta_i)'
            sum_alpha_sub_betas = np.zeros(D.shape)
            for i in range(self.n_respondents):
                alpha_sub_beta = alpha - beta[:, i]

                # logging.debug(alpha_sub_beta)
                # logging.debug(np.outer(alpha_sub_beta, alpha_sub_beta.T))

                sum_alpha_sub_betas += np.outer(alpha_sub_beta, alpha_sub_beta.T)

                # logging.debug(sum_alpha_sub_betas)
            # logging.debug("\n\n{}".format(sum_alpha_sub_betas))

            H = D + sum_alpha_sub_betas

            # logging.debug("sum_alpha_sub_betas: {}".format(sum_alpha_sub_betas))
            # logging.debug("H: {}".format(H))

            # H^-1 = TT'
            big_T = np.linalg.cholesky(np.linalg.inv(H))
            big_T = np.triu(big_T.T, 1) + big_T

            # logging.debug("big T: {}".format(big_T))

            # {u_1, ... , u_N | u_i ~ Normal(0, 1)} 
            # S = sum_N (T u_i)(T u_i)'
            u_vectors = SPECIAL_STATE.multivariate_normal(
                    np.zeros(self.n_parameters-self.n_attributes+1), 
                    np.eye(self.n_parameters-self.n_attributes+1), 
                    size=big_N)
            S = sum([np.outer(np.dot(big_T, u_vectors[i, :]), 
                                np.dot(big_T, u_vectors[i, :])) 
                        for i in range(big_N)])

            # logging.debug("S: {}".format(S))

            # D = S^-1
            D = np.linalg.inv(S)
            
            # logging.debug("D: {}".format(D))

            # Estimate Beta
            beta_old = beta
            differences = SPECIAL_STATE.multivariate_normal(
                    np.zeros(self.n_parameters-self.n_attributes+1), 
                    np.multiply(D, proportionality_factor), 
                    size=self.n_respondents)
            beta_new = beta + differences.T
            updated_beta = np.empty(beta.shape)

            # logging.debug("proportionality factor: {}".format(proportionality_factor))
            # logging.debug("Differences: {}".format(differences))

            # Count of respondents who accepted the new betas
            accepted = 0
            bad_ones = 0

            for respondent in range(self.n_respondents):
                try:
                    # logging.debug(self._calcProbabilityForRespondent(beta_new[:, respondent], respondent))
                    # logging.debug(self._calcDensity(alpha, D, beta_new[:, respondent]))
                    # logging.debug(self._calcProbabilityForRespondent(beta_old[:, respondent], respondent))
                    # logging.debug(self._calcDensity(alpha, D, beta_old[:, respondent]))

                    prob_new = self._calcProbabilityForRespondent(beta_new[:, respondent], respondent)
                    den_new = self._calcDensity(alpha, D, beta_new[:, respondent])
                    prob_old = self._calcProbabilityForRespondent(beta_old[:, respondent], respondent)
                    den_old = self._calcDensity(alpha, D, beta_old[:, respondent])
                except Exception as e:
                    raise ValueError("Error when calculating r \n\t{}".format(e))

                # Accept new beta if denominator is zero
                if (den_old == 0 or prob_old == 0):
                    updated_beta[:, respondent] = beta_old[:, respondent]
                    bad_ones += 1
                else: 
                    prob_ratio = prob_new / prob_old
                    den_ratio = den_new / den_old
                    r = prob_ratio * den_ratio

                    # logging.debug("r: {}".format(r))

                    # Accept new beta if r is greater than 1
                    if (r >= 1):
                        accepted += 1
                        updated_beta[:, respondent] = beta_new[:, respondent]
                    else:
                        # Accept new beta with probablility r
                        rand_val = SPECIAL_STATE.uniform()
                        if (rand_val <= r):
                            accepted += 1
                            updated_beta[:, respondent] = \
                                    beta_new[:, respondent]
                        else:
                            updated_beta[:, respondent] = \
                                    beta_old[:, respondent]

            beta = updated_beta

            logging.debug("beta: {}".format(beta))

            # Calculate acceptance rate and keep it around the ACCEPTANCE_RATE
            # Increase by 10% if greater than acceptance
            # Decrease by 10% if less than acceptance
            accepted_rate = 1.0 * accepted / self.n_respondents
            change_size = 0.1 * proportionality_factor
            if (accepted_rate > acceptance_rate):
                proportionality_factor += change_size
            if (accepted_rate < acceptance_rate):
                proportionality_factor -= change_size

            # logging.debug("*****************\nbad ones: {}\n***************".format(bad_ones))
            # logging.debug("accepted rate: {}".format(accepted_rate))

            # Store the alpha, beta, and avg_beta values used in this iteration
            self.alphas[iter_num, :] = alpha
            self.betas[iter_num, :, :] = beta
            self.avg_betas[iter_num, :] = avg_beta


    def _calcProbabilityForRespondent(self, beta, respondent):
        """Calculates the probabilites based off of the current beta 
            (partworths)

        Parameters
        ----------
        beta : Numpy array
            Numpy array containing the partworths of included parameters for
            each respondent and one for the normalizing term
            Dimensions are (n_parameters-n_attributes+1, n_respondents)
        respondent : int
            Column number of the respondent we are calculating the probability

        Returns
        -------

        
        """
        participated = np.where(self.selections.iloc[:, respondent] >= 0)
        resp_choices = lemonade._choices[participated, :, :]
        resp_options = lemonade._options[:, participated, :]
        n_participated = participated[0].shape[0]

        prod = 1
        for comp in range(n_participated):
            try:
                # Not sure why there is an extra dimension before the 
                # comparisons one but there is
                numer = np.exp(beta @ resp_choices[0, comp, respondent, :])
                denom = sum([np.exp(beta @ resp_options[i, 0, comp, :]) 
                                for i in range(self.n_options)])
                prod *= numer / denom
            except FloatingPointError:
                return 0
            except Exception as e:
                raise ValueError("Error in calculating probability \n\t{}".format(e))
                #logging.debug("comp: {}".format(comp))
                #logging.debug("Beta: {}".format(beta))
                #logging.debug("choice: {}".format(self.choices[comp, respondent, :]))
                #logging.debug("denom: {}".format(denom))
            
        # logging.debug("Probability: {}".format(prod))

        return prod


    def _calcDensity(self, alpha, D, beta):
        """Calculates the relative density

        Parameters
        ----------
        alpha : Numpy array
            Numpy array containing the mean partworth of included parameters
            and one for the normalizing term
            Dimensions are (n_parameters-n_attributes+1)
        beta : Numpy array
            Numpy array containing the partworths of included parameters for
            each respondent and one for the normalizing term
            Dimensions are (n_parameters-n_attributes+1, n_respondents)
        D : Numpy array
            Numpy array containing the covarience matrix for the betas
            Dimensions are (n_parameters-n_attributes x n_parameters-n_attributes)

        Returns
        -------

        
        """
        beta_sub_alpha = beta - alpha
        inv_D = np.linalg.inv(D)

        try:
            result = np.exp(-0.5 * beta_sub_alpha.T @ inv_D @ beta_sub_alpha)
        except FloatingPointError:
            result = 0

        # logging.debug("Density: {}".format(result))

        return result


    def calcRelativeImportance(self):
        """Calculates the relative importance of the results
        
        Yields
        ------
        self.relativeImportance : Numpy array
            Vector of the relative importances of each attribute

        Parameters
        ----------

        Returns
        -------
        Returns a numpy array containing the relative importance of each 
        attribute. Dimensions are (n_attributes)

        """
        ranges = np.zeros((self.n_attributes, self.n_respondents))

        for attr in range(self.n_attributes):
            for resp in range(self.n_respondents):
                values = list()
                values.append(self.betas[self.iters-1, 0, resp])
                
                ind = int(self._removed[attr])
                for i in range(int(self.n_attr_parameters[attr]) - 1):
                    values.append(self.betas[self.iters-1, ind-attr+1+i, resp])

                ranges[attr, resp] = max(values) - min(values)

        self.relative_importance = np.sum(ranges, axis=1) / np.sum(ranges)

        return self.relative_importance


if __name__ == "__main__":
    # Import in test data
    testfolder = "immigration_cjoint"
    sep = "\t"
    profiles = pd.read_csv("data/{}/profiles.tsv".format(testfolder), sep=sep).set_index('Profile')
    comparisons = pd.read_csv("data/{}/comparisons.tsv".format(testfolder), sep=sep).set_index('Comparisons')
    selections = pd.read_csv("data/{}/selections.tsv".format(testfolder), sep=sep).set_index("Comparisons") - 1
    lemonade = Cbc_hb(profiles, comparisons, selections)

    # Run 10 iterations 20,000 times each to determine convergence
    repetitions = 10
    iters = 1000
    repers = [i for i in range(repetitions)]
    alphas = np.zeros((repetitions, iters, lemonade.n_parameters - lemonade.n_attributes))
    avg_betas = np.zeros((repetitions, iters, lemonade.n_parameters - lemonade.n_attributes))
    betas = np.zeros((repetitions, iters, lemonade.n_parameters - lemonade.n_attributes, lemonade.n_respondents))
    probabilities = np.zeros((repetitions, lemonade.n_parameters))
    relative_importance = np.zeros((repetitions, lemonade.n_attributes))

    # Run multiple iterations
    if False:
        while repers:
            i = repers[0]
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            try:
                print("iter {} started".format(i))
                lemonade.model(iters=iters)
                alphas[i, :, :] = lemonade.alphas
                avg_betas[i, :, :] = lemonade.avg_betas
                betas[i, :, :, :] = lemonade.betas
                relative_importance[i, :] = lemonade.calcRelativeImportance()
                print("iter {} passed".format(i))
                repers.remove(i)
            except Exception as e:
                print("iter {} failed\n******\n{}\n******".format(i, e))
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        np.save('tmp/alphas.npy', alphas)
        np.save('tmp/avg_betas.npy', avg_betas)
        np.save('tmp/betas.npy', betas)
        np.save('tmp/relative_importance.npy', relative_importance)

    # Run one iteration
    if True:
        failed = True
        while(failed):
            try:
                lemonade.model(iters=iters)
                lemonade.calcRelativeImportance()
                failed = False
            except Exception as e:
                failed = True
                break

        np.save('tmp/alphas.npy', lemonade.alphas)
        np.save('tmp/avg_betas.npy', lemonade.avg_betas)
        np.save('tmp/betas.npy', lemonade.betas)
        np.save('tmp/relative_importance.npy', lemonade.relative_importance)

