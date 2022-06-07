# The functions within this python file are for generating predictions of 
#   protodeboronation using your own DFT and pKa/pKaH data. Parameters for 
#   the linear sklearn models are loaded from a pickle file.
# Constructor takes the DFT calculated energy for the active mechanistic
#   pathways, all other mechanistic pathways should be set to None
#       Example: k1 = 0.0152523, k2 = -0.03231252, k2Ar = None, (...)
import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle


class BoronicAcid:
    def __init__(self, k1, k2, k2Ar, k2cat, k3, k4, k5, pKa, pKaH=None ):
        self.k1 = self._mechanism_rate(k1, 'k1')
        self.k2 = self._mechanism_rate(k2, 'k2')
        self.k2Ar = self._mechanism_rate(k2Ar, 'k2Ar')
        self.k2cat = self._mechanism_rate(k2cat, 'k2cat') 
        self.k3 = self._mechanism_rate(k3, 'k3')
        self.k4 = self._mechanism_rate(k4, 'k4')
        self.k5 = self._mechanism_rate(k5, 'k5')
        self.pKa = pKa
        self.pKaH = pKaH

    def _mechanism_rate(self, delta_E, mechanism):
        """
        Predict the max rate for a particular mechanism, given delta_E of that mechanism
        """
        # if delta_E is none, return nothing, else calculate the mechanism rate
        if delta_E != None:
            # Load model
            filename = f"models/{mechanism}_model.sav"
            loaded_model = pickle.load(open(filename, 'rb'))
            

            # predicted max rate of this mechanism
            predicted_max_rate = loaded_model.predict(
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

        if self.k1 != None:
            rate += 10 ** (self.k1 - 1 * pH)

        if self.k2 != None:
            if 0 <= pH and pH < self.pKa:
                b = (
                    self.k2 - 0.75 * self.pKa
                )  # -0.75 was arrived at by inspection of the data, see the paper for more details
                rate += 10 ** (0.75 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                rate += 10 ** (self.k2)

        if self.k2Ar != None:
            if 0 <= pH and pH < self.pKa:
                # y = ax+b => b = y - ax => b = k_max - 0.75*pKa
                b = self.k2Ar - 0.75 * self.pKa
                rate += 10 ** (0.75 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                rate += 10 ** (self.k2Ar)

        if self.k2cat != None:
            if 0 <= pH and pH < self.pKa:
                # y = ax+b => b = y - ax => b = k_max - 2 * pKa
                b = self.k2cat - 2 * self.pKa
                rate += 10 ** (2 * pH + b)
            elif self.pKa <= pH and pH <= 14:
                # y = ax+b => b = y - ax => b = k_max - (-2) * pKa
                b = self.k2cat + 2 * self.pKa
                rate += 10 ** (-2 * pH + b)

        if self.k3 != None:
            # y=ax+b => b = y-ax => b = k_max - 2 * max_pH
            b = self.k3 - 2 * 14
            rate += 10 ** (2 * pH + b)

        if self.k4 != None:
            if 0 <= pH and pH < self.pKaH:
                # y = ax+b => b = y - ax => b = k_max - 0.75*pKaH
                b = self.k4 - 0.75 * self.pKaH
                rate += 10 ** (0.75 * pH + b)
            elif self.pKaH <= pH and pH <= self.pKa:
                rate += 10 ** (self.k4)
            elif self.pKa < pH and pH <= 14:
                b = self.k4 + 0.75 * self.pKa
                rate += 10 ** (-0.75 * pH + b)

        if self.k5 != None:
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
