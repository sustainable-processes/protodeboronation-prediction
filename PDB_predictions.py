# The functions within this python file are for generating predictions for the rate 
#   of protodeboronation, either using your own DFT and pKa/pKaH data 
#   (BoronicAcid parent class), or using the data found in the csv files of this 
#   directory (allowing you to reproduce the results found in the paper) 
#   (KnownBoronicAcid child class). 
# Parameters for the linear sklearn models are loaded from a pickle file when
#   BoronicAcid is used, and trained anew/in real time from the csv data when 
#   KnownBoronicAcid is used.
#
# The BoronicAcid constructor takes the DFT calculated energy for the 
#   active mechanistic pathways, all other mechanistic pathways should 
#   be set to None.
#       Example: k1 = 0.0152523, k2 = -0.03231252, k2Ar = None, (...)
#
# The KnownBoronicAcid constructor simply takes the molecular ID of a molecule
#   (which can be found in the paper), and automatically extracts the relevant 
#   data from the csv files.
#

import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle


class BoronicAcid:
    def __init__(self, k1, k2, k2Ar, k2cat, k3, k4, k5, pKa, pKaH, smiles=None):
        self.k1 = self._mechanism_rate(k1, 'k1')
        self.k2 = self._mechanism_rate(k2, 'k2')
        self.k2Ar = self._mechanism_rate(k2Ar, 'k2Ar')
        self.k2cat = self._mechanism_rate(k2cat, 'k2cat') 
        self.k3 = self._mechanism_rate(k3, 'k3')
        self.k4 = self._mechanism_rate(k4, 'k4')
        self.k5 = self._mechanism_rate(k5, 'k5')
        self.pKa = pKa
        self.pKaH = pKaH #set to None if there is no pKaH value
    
    def __repr__(self) -> str:
        return "BoronicAcid({}, {}, {}, {}, {}, {}, {}, {}, {})".format(self.k1, self.k2, self.k2Ar, self.k2cat, self.k3, self.k4, self.k5, self.pKa, self.pKaH)

    # def __str__(self) -> str:
    #     pass

    def _get_model(self, mechanism):
        filename = f"models/{mechanism}_model.sav"
        loaded_model = pickle.load(open(filename, 'rb'))
        return loaded_model

    def _mechanism_rate(self, delta_E, mechanism):
        """
        Predict the max rate for a particular mechanism, given delta_E of that mechanism
        """
        # if delta_E is none, return nothing, else calculate the mechanism rate
        if delta_E is not None:

            model = self._get_model(mechanism)
            

            # predicted max rate of this mechanism
            predicted_max_rate = model.predict(
                np.array(delta_E).reshape(1, -1)
            )
            predicted_max_rate = predicted_max_rate[0]

            # rate cannot be below -8 or above 2
            if predicted_max_rate < -8:
                predicted_max_rate = -8
            elif predicted_max_rate > 2:
                predicted_max_rate = 2

            return predicted_max_rate

        else:
            return None

    def total_rate_point_predictor(self, pH: float):
        """
        The total (observed) rate of protodeboronation at a particular pH
            is the sum of rates for all active mechanistic pathways
        The shape of each mechanistic rate curve is explained in the paper
        This function returns the predicted total rate for self.molecule_ID at pH
        """
        # initialise the rate
        rate = 0.0

        if self.k1 is not None:
            rate += 10 ** (self.k1 - 1 * pH)

        if self.k2 is not None:
            if 0 <= pH and pH < self.pKa:
                b = (
                    self.k2 - 0.75 * self.pKa
                )  # -0.75 was arrived at by inspection of the data, see the paper for more details
                rate += 10 ** (0.75 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                rate += 10 ** (self.k2)

        if self.k2Ar is not None:
            if 0 <= pH and pH < self.pKa:
                # y = ax+b => b = y - ax => b = k_max - 0.75*pKa
                b = self.k2Ar - 0.75 * self.pKa
                rate += 10 ** (0.75 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                rate += 10 ** (self.k2Ar)

        if self.k2cat is not None:
            if 0 <= pH and pH < self.pKa:
                # y = ax+b => b = y - ax => b = k_max - 2 * pKa
                b = self.k2cat - 2 * self.pKa
                rate += 10 ** (2 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                # y = ax+b => b = y - ax => b = k_max - (-2) * pKa
                b = self.k2cat + 2 * self.pKa
                rate += 10 ** (-2 * pH + b)

        if self.k3 is not None:
            # y=ax+b => b = y-ax => b = k_max - 2 * max_pH
            b = self.k3 - 2 * 14
            rate += 10 ** (2 * pH + b)

        if self.k4 is not None:
            if 0 <= pH and pH < self.pKaH:
                # y = ax+b => b = y - ax => b = k_max - 0.75*pKaH
                b = self.k4 - 0.75 * self.pKaH
                rate += 10 ** (0.75 * pH + b)
            elif self.pKaH <= pH and pH <= self.pKa:
                rate += 10 ** (self.k4)
            elif self.pKa < pH and pH <= 14:
                b = self.k4 + 0.75 * self.pKa
                rate += 10 ** (-0.75 * pH + b)

        if self.k5 is not None:
            if 0 <= pH and pH < self.pKaH:
                rate += 10 ** (self.k5)
            elif self.pKaH <= pH and pH <= 14:
                # y = ax+b => b = y - ax => b = k_max - 0.75*pKaH
                b = self.k5 + 0.75 * self.pKaH
                rate += 10 ** (-0.75 * pH + b)
        if rate == 0:
            k_obs = -9
        else:
            k_obs = np.log10(rate)
        if k_obs < -9:
            k_obs = -9
        if k_obs > 2:
            k_obs = 2

        return k_obs

    def plot_result(self,title):
        """
        Plot the predicted (solid blue) and measured (green dots) log(k) vs pH

        Measured data is only available for the 50 Cox molecules, so the measured
            rate will obviously only be plotted when choosing one of these 50 molecules
        """

        fig_size = (10, 5)
        f = plt.figure(figsize=fig_size)

        X_pred = []
        y_pred = []
        i = np.linspace(0, 14, 140)
        for pH in i:
            X_pred += [pH]
            y_pred += [self.total_rate_point_predictor(pH)]

        plt.plot(X_pred, y_pred, "b-", label="Predicted rate")
        plt.legend(prop={"size":12})
        plt.ylim([-10,2])

        plt.title(title, fontdict={"fontsize": 18})
        plt.legend(prop={"size": 12})
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylim([-10, 2.2])
        plt.xlabel("pH", fontdict={"fontsize": 14})
        plt.ylabel("log(k)", fontdict={"fontsize": 14})
        plt.tight_layout()

        return f
    
    # Source:
    # https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Nuclear_Chemistry/Nuclear_Kinetics/Half-Lives_and_Radioactive_Decay_Kinetics
    def halflife(self,pH):
        log_k = self.total_rate_point_predictor(pH)
        k = 10**log_k #rate 
        t = math.log(2)/k #half life
        
        t_mseconds = t/1000
        t_seconds = t
        t_minutes = t/60
        t_hours = t/(60*60)
        t_days = t/(60*60*24)
        t_months = t_days/30
        
        #choose an appropriate unit:
        if t_seconds < 1:
            return round(t_mseconds,5), 'milliseconds'
        if t_seconds < 120:
            return round(t_seconds,2), 'seconds'
        elif t_minutes < 120:
            return round(t_minutes,2), 'minutes'
        elif t_hours < 48:
            return round(t_hours,2), 'hours'
        elif t_days < 60:
            return round(t_days,2), 'days'
        else:
            return round(t_months,2), 'months'

class KnownBoronicAcid(BoronicAcid):
    #This class allows us to reproduce the results shown in the paper by pulling
    #   data from the csv file, using just molecule ID
    def __init__(self, molecule_ID):
        self.molecule_ID = molecule_ID
        
        # set all the mechanism rates
        self.k1 = self._mechanism_rate(self._extract_delta_E("k1"), "k1")
        self.k2 = self._mechanism_rate(self._extract_delta_E("k2"), "k2")
        self.k2Ar = self._mechanism_rate(self._extract_delta_E("k2Ar"), "k2Ar")
        self.k2cat = self._mechanism_rate(self._extract_delta_E("k2cat"), "k2cat")
        self.k3 = self._mechanism_rate(self._extract_delta_E("k3"), "k3")
        self.k4 = self._mechanism_rate(self._extract_delta_E("k4"), "k4")
        self.k5 = self._mechanism_rate(self._extract_delta_E("k5"), "k5")

        # set the pKa and pKaH (if applicable)
        # pKa and pKaH info:
        df = pd.read_csv("data/Cox-molecules/Cox-molecules-overview.csv")

        # set pKa and pKaH values
        if self.molecule_ID < 100:  # It's a Cox molecule
            # remove all rate measurements that aren't about the query molecule
            df = df[df.molecule_number == self.molecule_ID]
            self.pKa = list(df["pKa"])[0]
            self.pKaH = list(df["pKaH"])[0]
            self.smiles = list(df["smiles"])[0]
        elif self.molecule_ID >= 100 and self.molecule_ID <= 150:  # It's a novel molecule
            # simply use the average pKa and pKaH of all the Cox molecules
            pKa_series = df["pKa"].dropna()
            self.pKa = sum(list(pKa_series)) / len(list(pKa_series))
            pKaH_series = df["pKaH"].dropna()
            self.pKaH = sum(list(pKaH_series)) / len(list(pKaH_series))

            novel_df = pd.read_csv("data/Novel-molecules/Novel-molecules-data.csv")
            self.smiles = list(novel_df["smiles"])[0]
        
        if self.molecule_ID not in list(df["molecule_number"]) and self.molecule_ID not in range(100,151):
            raise IndexError("Molecule_ID not valid")

    class IndexError(Exception):
        pass

    
    def __repr__(self) -> str:
        return "KnownBoronicAcid({})".format(self.molecule_ID)
    
    def __str__(self) -> str:
        return self.smiles

    def _extract_delta_E(self, mechanism: str):
        """
        Each mechanism follow the scheme transition state -> products
        The energy of each molecule/molecular structure was carried out using DFT
        Taking the difference in energy between the products and the transition state yields delta_E
        This energy difference was saved in a csv file
        Each boronic acid can degrade along one or more mechanistic pathways

        This function extracts delta_E for a particular mechanism (e.g. 'k1')
            of self.molecule_ID and returns delta_E: float

        """
        # Only these 7 mechanisms are defined for boronic acids of type -B(O)O
        assert mechanism in ["k1", "k2", "k2Ar", "k2cat", "k3", "k4", "k5"]

        if self.molecule_ID < 100:  # It's a Cox molecule
            df = pd.read_csv(f"data/Cox-molecules/{mechanism}.csv")
        elif self.molecule_ID >= 100:  # It's a Novel molecule
            df = pd.read_csv(f"data/Novel-molecules/{mechanism}.csv")
        delta_E = df.loc[df["molecule_number"] == self.molecule_ID, "delta_E"]
        delta_E = list(delta_E)

        try:
            delta_E = delta_E[0]
            return delta_E
        except IndexError:
            return None

    def _get_model(self, mechanism):
        # Train the model
        df = pd.read_csv(f"data/Cox-molecules/{mechanism}.csv")
        # remove molecule_ID from the dataset
        # i.e. ensure that data about the test molecule isn't used for training
        training_df = df[df["molecule_number"] != self.molecule_ID]

        # also remove any rows containing null values
        training_df = training_df[training_df.delta_E.notnull()]
        training_df = training_df[training_df.k_obs_max.notnull()]

        # extract the transition state energy and rate data and transform to np.array
        delta_E_training = np.array(training_df["delta_E"]).reshape((-1, 1))
        k_obs_training = np.array(training_df["k_obs_max"])

        # train linear regression model
        model = LinearRegression().fit(
            delta_E_training, k_obs_training
        )
        return model
    
    def plot_result(self):

        title = self.molecule_ID

        f = super().plot_result(title)

        if self.molecule_ID < 100:
            # remove all rate measurements that aren't about the query molecule
            df = pd.read_csv(f"data/Cox-molecules/Cox-molecules-data.csv")
            query_molecule_df = df[df.molecule_number == self.molecule_ID]

            X = query_molecule_df["pH"]
            y = query_molecule_df["log(k_obs)"]
            plt.plot(X, y, "go", label="Measured rate")
        
        return f

        
        