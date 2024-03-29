

# Fracture energy calculation

import numpy as np


class FractureEnergy:

    def __init__(self, G0, Lcnt, Dcnt, SigmaUlt, TauInt, Ac, mu,
                 Ecnt, ThetaMin, ThetaMax, fc, p=0.5, q=0.5,
                 NInter=100, Ntheta=400):
        #
        # ------------------------------------------------------ #
        # -----------------Fracture Energy---------------------- #
        # ------------------------------------------------------ #
        #
        # Fracture energy of the pristine
        self.G0 = G0
        # Cnt volume fraction
        self.fc = fc
        # Carbon nanotube length
        self.Lcnt = Lcnt
        # Carbon nanotube diameter
        self.Dcnt = Dcnt
        # Ultimate stress
        self.SigmaUlt = SigmaUlt
        # Friction shear stress
        self.TauInt = TauInt
        # A coeff
        self.Ac = Ac
        # Area
        self.A = np.power(Dcnt, 2) / 4
        # Snubbing coeff
        self.mu = mu
        # Cnt Young's modulus
        self.Ecnt = Ecnt
        # Theta min integration
        self.ThetaMin = ThetaMin
        # Theta min integration
        self.ThetaMax = ThetaMax
        # Discretization of the integration
        # of g(theta)
        self.NInter = NInter
        # Discretization of the integration
        # of theta
        self.Ntheta = Ntheta
        # p
        self.p = p
        # q
        self.q = q

    def Lctheta(self, theta):
        #
        # ------------------------------------------------------ #
        # -----------------Critical Length---------------------- #
        # ------------------------------------------------------ #
        #
        num = self.Dcnt * (self.SigmaUlt * (1 - self.Ac * np.tan(theta)))

        den = self.TauInt * np.exp(self.mu * theta) * 4

        return num / den

    def gtheta(self, theta):
        #
        # ------------------------------------------------------ #
        # -----------------g(theta) function-------------------- #
        # ------------------------------------------------------ #
        #
        Nint = self.NInter

        def num(theta): return (np.sin(theta)**(2 * self.p - 1.)) * \
            (np.cos(theta)**(2 * self.q - 1.))

        x = np.linspace(self.ThetaMin, self.ThetaMax, Nint)

        int1 = np.zeros([Nint])

        for i in range(len(x)):

            int1[i] = num(x[i])

        den = np.trapz(int1, x)

        return num(theta) / den

    def Integranddl(self, theta):
        #
        # ------------------------------------------------------ #
        # --------------------Integrand dl---------------------- #
        # ------------------------------------------------------ #
        #
        LcTheta = self.Lctheta(theta)

        if LcTheta <= self.Lcnt:

            value = np.power(LcTheta, 3) * self.TauInt * np.pi * \
                self.Dcnt * np.exp(self.mu * theta) / 48 + \
                (self.Lcnt - LcTheta) * np.pi * np.power(self.Dcnt, 2) * \
                np.power(self.SigmaUlt, 2) * self.Lcnt / (16 * self.Ecnt)

        else:

            value = np.power(self.Lcnt, 3) * self.TauInt * np.pi * \
                self.Dcnt * np.exp(self.mu * theta) / 48

        if value < 0:

            return 0

        else:

            return value

    def EnergyReleaseRate(self):

        Nint = self.Ntheta

        ThetaSerie = np.linspace(self.ThetaMin, self.ThetaMax, Nint)

        Integrand = np.zeros([Nint])

        for i in range(Nint):

            theta = ThetaSerie[i]

            factor = (2 * self.fc) / (self.A * self.Lcnt)

            Integrand[i] = factor * self. gtheta(theta) * np.cos(theta) * \
                self. Integranddl(theta)

        return np.trapz(Integrand, ThetaSerie) + self.G0


if __name__ == "__main__":

    Lcnt = 3.20995854347668e-6  # Length [microns]

    Dcnt = 10.3538294770233e-9  # Diameter [nm]

    ThetaMin = 0  # [rad]

    ThetaMax = np.pi / 2  # [rad]

    LMin = 0

    LMax = Lcnt / 2  # Maximum Cnt length to integrate

    mu = 0  #
    
    Ac = 0.083  # Orientation limit angle
    
    G0 = 200  # [J/m2]
    
    SigmaUlt = 35e9  # CNT ultimate strength [Pa]
    
    Ecnt = 700e9  # Youngs modulus MPa
    
    nu_cnt = 0.3  # Poissons ratio
    
    TauInt = 47e6  # CNT interfatial strength [Pa] (epoxy-CNT)
    
    fc = 0.01
    
    GM = FractureEnergy(G0, Lcnt, Dcnt, SigmaUlt, TauInt, Ac, mu,
                        Ecnt, ThetaMin, ThetaMax, fc, p=0.5, q=0.5,
                        NInter=100, Ntheta=1000)

    FE = GM.EnergyReleaseRate()

    print(FE)
