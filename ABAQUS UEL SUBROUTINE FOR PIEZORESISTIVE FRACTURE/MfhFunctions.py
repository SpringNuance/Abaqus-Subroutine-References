# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 23:12:24 2023

@author: Leonel Quinteros
"""

import numpy as np


class MFH:

    def __init__(self, Lcnt, Dcnt, EMatrix, NuMatrix, ECnt,
                 NuCnt, EInt, NuInt, Intert, fc):
        #
        # ------------------------------------------------------ #
        # ----------------Elastic properties-------------------- #
        # ------------------------------------------------------ #
        #
        self.Lcnt = Lcnt

        self.Dcnt = Dcnt

        self.EMatrix = EMatrix

        self.NuMatrix = NuMatrix

        self.ECnt = ECnt

        self.NuCnt = NuCnt

        self.EInt = EInt

        self.NuInt = NuInt

        self.Intert = Intert

        self.fc = fc

    def doubleintegral_simpson(self, function_phi_theta, phi, theta):
        #
        # ------------------------------------------------------ #
        # ---------------DoubleInteral simpson------------------ #
        # ------------------------------------------------------ #
        #
        number_phi = phi.shape[0]

        number_theta = theta.shape[0]

        lowlimit_phi = phi[0]

        highlimit_phi = phi[-1]

        lowlimit_theta = theta[0]

        highlimit_theta = theta[-1]

        d_phi = (highlimit_phi - lowlimit_phi) / (number_phi - 1)

        d_theta = (highlimit_theta - lowlimit_theta) / (number_theta - 1)

        # Integral over theta

        Simpsonmatrix_theta = np.zeros([number_theta])
        Simpsonmatrix_theta[0] = 1
        Simpsonmatrix_theta[2:-2:2] = 2
        Simpsonmatrix_theta[1:-1:2] = 4
        Simpsonmatrix_theta[-1] = 1
        theta_int = np.mat(function_phi_theta) * np.mat(Simpsonmatrix_theta).T

        # Integral over phi
        Simpsonmatrix_phi = np.zeros([number_phi])
        Simpsonmatrix_phi[0] = 1
        Simpsonmatrix_phi[2:-2:2] = 2
        Simpsonmatrix_phi[1:-1:2] = 4
        Simpsonmatrix_phi[-1] = 1
        b = np.mat(Simpsonmatrix_phi) * theta_int
        phi_theta_int_int = (d_theta * d_phi) * b / 9

        return phi_theta_int_int[0, 0]

    def computeEngineeringConstantsSqrt2(self, C):

        E = np.zeros([6])
        nu = np.zeros([6])

        # if (rank(C) < 6)
        #  return;
        # end

        # Compute compliance matrix:
        # Inverse of a symmetric matrix is symmetric
        S = np.mat(np.linalg.inv(C))

        for i in range(6):
            E[i] = 1.0 / S[i, i]

        # Define nu: [nu_23 nu_13 nu_12 nu_32 nu_31 nu_21]

        # Better accuracy for larger poisson ratio:
        nuscale = -(np.diag(E[0:3]) * np.transpose(S[0:3, 0:3]))

        nu[0] = nuscale[1, 2]
        nu[1] = nuscale[0, 2]
        nu[2] = nuscale[0, 1]
        nu[3] = nuscale[2, 1]
        nu[4] = nuscale[2, 0]
        nu[5] = nuscale[1, 0]
        return E, nu

    def int_mat(self, M_to_int, phi, theta):
        # Compute integration of matrix by the Simpson Method
        #
        # INPUT:
        # M_to_int -> Matrix to be integrated
        # phi -> Array of values for rotation angle phi
        # theta -> Array of values for rotation angle theta
        #
        # OUTPUT:
        # Mav -> Integrated Matrix

        # Author: E. García-Macías
        # ------------------------------------------------------------------------
        length_theta = theta.shape[0]
        length_phi = phi.shape[0]
        d1, d2, d3, d4 = M_to_int.shape
        Mav = np.zeros([d1, d2])
        for i in range(d1):
            for j in range(d2):
                Mav[i, j] = self.doubleintegral_simpson(np.reshape(
                    M_to_int[i, j, :, :], [length_phi, length_theta]), phi, theta)
        return Mav

    def rotmatrix(self, Beta, Alpha, Psi, var):

        # eul=[Psi Alpha Beta];
        Beta = Beta * 180 / np.pi
        Alpha = Alpha * 180 / np.pi
        Psi = Psi * 180 / np.pi

        if var == []:
            orden = [3, 2, 1]
        else:
            orden = var

        R = np.zeros([3, 3, 3])

        R[0, :, :] = [[1, 0, 0], [0, np.cos(np.radians(Beta)), np.sin(np.radians(Beta))], [
            0, -np.sin(np.radians(Beta)), np.cos(np.radians(Beta))]]
        R[1, :, :] = [[np.cos(np.radians(Alpha)), 0, -np.sin(np.radians(Alpha))],
                      [0, 1, 0], [np.sin(np.radians(Alpha)), 0, np.cos(np.radians(Alpha))]]
        R[2, :, :] = [[np.cos(np.radians(Psi)), np.sin(np.radians(
            Psi)), 0], [-np.sin(np.radians(Psi)), np.cos(np.radians(Psi)), 0], [0, 0, 1]]

        A1 = np.mat(np.squeeze(R[orden[2], :, :]))
        A2 = np.mat(np.squeeze(R[orden[1], :, :]))
        A3 = np.mat(np.squeeze(R[orden[0], :, :]))

        w = (A1 * A2) * A3
        # w = np.matmul(np.matmul(np.squeeze(R[orden[2],:,:]),np.squeeze(R[orden[1],:,:])),np.squeeze(R[orden[0],:,:]))

        # eul2rotm(eul,'ZYX')

        K1 = np.array([[w[0, 0]**2, w[0, 1]**2, w[0, 2]**2], [w[1, 0]**2,
                      w[1, 1]**2, w[1, 2]**2], [w[2, 0]**2, w[2, 1]**2, w[2, 2]**2]])
        K2 = np.array([[w[0, 1] *
                        w[0, 2], w[0, 2] *
                        w[0, 0], w[0, 0] *
                        w[0, 1]], [w[1, 1] *
                       w[1, 2], w[1, 2] *
                       w[1, 0], w[1, 0] *
                       w[1, 1]], [w[2, 1] *
                       w[2, 2], w[2, 2] *
                       w[2, 0], w[2, 0] *
                       w[2, 1]]])
        K3 = np.array([[w[1, 0] *
                        w[2, 0], w[1, 1] *
                        w[2, 1], w[1, 2] *
                        w[2, 2]], [w[2, 0] *
                       w[0, 0], w[2, 1] *
                       w[0, 1], w[2, 2] *
                       w[0, 2]], [w[0, 0] *
                       w[1, 0], w[0, 1] *
                       w[1, 1], w[0, 2] *
                       w[1, 2]]])
        K4 = np.array([[w[1, 1] *
                        w[2, 2] +
                        w[1, 2] *
                        w[2, 1], w[1, 2] *
                        w[2, 0] +
                        w[1, 0] *
                        w[2, 2], w[1, 0] *
                        w[2, 1] +
                        w[1, 1] *
                        w[2, 0]], [w[2, 1] *
                       w[0, 2] +
                       w[2, 2] *
                       w[0, 1], w[2, 2] *
                       w[0, 0] +
                       w[2, 0] *
                       w[0, 2], w[2, 0] *
                       w[0, 1] +
                       w[2, 1] *
                       w[0, 0]], [w[0, 1] *
                       w[1, 2] +
                       w[0, 2] *
                       w[1, 1], w[0, 2] *
                       w[1, 0] +
                       w[0, 0] *
                       w[1, 2], w[0, 0] *
                       w[1, 1] +
                       w[0, 1] *
                       w[1, 0]]])

        gchange3D1 = np.mat(np.zeros([6, 6]))

        gchange3D1[0:3, 0:3] = K1

        gchange3D1[0:3, 3:6] = 2 * K2

        gchange3D1[3:6, 0:3] = K3

        gchange3D1[3:6, 3:6] = K4

        gchange3D2 = np.mat(np.zeros([6, 6]))

        gchange3D2[0:3, 0:3] = K1
        
        gchange3D2[0:3, 3:6] = K2
        
        gchange3D2[3:6, 0:3] = 2 * K3
        
        gchange3D2[3:6, 3:6] = K4

        return gchange3D1, gchange3D2

    def rottensor(self, Beta, Alpha, Psi, order, A, val):
        # W -> Matrix to facilitate transformation between matrix and
        # tensor notation

        W = np.zeros([6, 6])

        W[0, 0] = 1

        W[1, 1] = 1

        W[2, 2] = 1

        W[3, 3] = 2

        W[4, 4] = 2

        W[5, 5] = 2

        invW = np.zeros([6, 6])

        invW[0, 0] = 1

        invW[1, 1] = 1

        invW[2, 2] = 1

        invW[3, 3] = 0.5

        invW[4, 4] = 0.5

        invW[5, 5] = 0.5

        T1, T2 = self.rotmatrix(Beta, Alpha, Psi, order)

        T1T = np.mat(np.transpose(T1))

        T2T = np.mat(np.transpose(T2))

        if val == 1:  # Stiffness tensor

            A_rot = np.linalg.inv(T2T) * (A * T1T)

        elif val == 2:  # Compliance tensor

            A_rot = np.linalg.inv(T1T) * (A * T2T)

        elif val == 3:  # Eshelby tensor

            A_rot = T2 * (A * np.linalg.inv(T2))

        elif val == 4:  # B

            A_rot = T1 * (A * np.linalg.inv(T1))

        return A_rot

    def orientational_average(self, number_theta, number_phi,
                              theta_max, phi_max, A, k, val):

        # Compute orientational average of tensor A
        #
        # INPUT:
        # numer_theta -> number of steps of integration for the variable theta
        # number_phi -> number of steps of integration for the variable phi
        # theta_max -> maximum value for variable theta
        # phi_max -> maximum value for variable phi
        # A -> Tensor to be averaged
        # k -> Parameter of the ODF according to "Constitutive modeling of
        # nanotube–reinforced polymer composites" Odegard et al.
        #
        # OUTPUT:
        # Amed -> Orientational average of input tensor A
        #
        # ------------------------------------------------------------------------

        if k > 1000000:

            Amed = A

        else:

            theta = np.linspace(0, theta_max, number_theta + 1)

            phi = np.linspace(0, phi_max, number_phi + 1)

            # Orientation distribution function
            # -----------------------------------

            length_theta = number_theta + 1

            length_phi = number_phi + 1

            # caso=1 random
            # caso=2 aligned
            # caso=3 general

            caso = 3

            ODF = np.zeros([length_theta])

            if caso == 1:

                ODF = ODF + 1

            elif caso == 2:

                ODF[0] = 1

            else:
                # k=0; % Inf-> Aligned, 0-> Random

                ODF = np.exp(-k * theta**2)

            # INTEGRATION
            # -----------------------------------

            # Inititalize result matrixs

            A_phi_theta = np.zeros([6, 6, length_phi, length_theta])

            denom_phi_theta = np.zeros([1, 1, length_phi, length_theta])

            A_to_int = A_phi_theta

            for t in range(length_theta):

                theta_loop = theta[t]

                for s in range(length_phi):

                    phi_loop = phi[s]

                    A_phi_theta[:, :, s, t] = self.rottensor(
                        theta_loop, phi_loop, 0, [2, 0, 1], A, val)

                    A_phi_theta[:, :, s, t] = self.rottensor(
                        0, phi_loop, theta_loop, [0, 2, 1], A, val)

                denom_phi_theta[0, 0, :, t] = (ODF[t] * np.sin(theta_loop))

                A_to_int[:, :, :, t] = np.reshape(A_phi_theta[:, :, :, t], [
                                                  6, 6, len(phi)]) * ODF[t] * np.sin(theta_loop)

            # Perform double integration

            denom = self.int_mat(denom_phi_theta, phi, theta)

            Amed = self.int_mat(A_to_int, phi, theta) / (denom)

            return np.mat(Amed)

    def softinterphase(self, kappa, lambd, vp):

        if vp >= 1:

            vi = 0

        else:

            if kappa > 1:

                eta = 1 / kappa

            else:

                eta = kappa

        phi = np.arccos(eta)

        if kappa < 1:

            nkappa = ((2 * kappa**(2. / 3.)) * np.sin(phi)) / \
                (np.sin(phi) + (kappa**2) * np.arctanh(np.sin(phi)))

        elif kappa == 1:

            nkappa = 1

        else:

            nkappa = ((2 * kappa**(2. / 3.)) * np.tan(phi)) / \
                (np.tan(phi) + (kappa**2) * phi)
        m = 0

        vi = (1 - vp) * (1 - np.exp(-(6 * vp / (1 - vp)) *
                                    (lambd / nkappa + (2 + 3 * vp / ((nkappa**2) * (1 - vp))) *
                                    (lambd**2) + (4. / 3.) * (1 + 3 * vp / (nkappa * ((1 - vp)**2)) +
                                                              m * (vp**2) / ((nkappa**3) * ((1 - vp)**2))) * lambd**3)))

        return vi

    def ellipsoidalinter_random2(self, Cp, Ci, Cm, Si, Sp, vp, vi):

        if vi == 0:

            number_theta = 20

            number_phi = 20

            theta_max = np.pi / 2

            phi_max = 2 * np.pi

            if vp == 0:

                Ceff = Cm

            else:

                Id = np.eye(6)

                vm = 1 - vp - vi

                psiP = Cm * np.linalg.inv(Cp - Cm)

                psiI = Cm * np.linalg.inv(Ci - Cm)

                dS = (Si - Sp)

                factor1 = (Sp + psiP - (vp / vi) * dS)

                factor2 = (Sp + psiI - (vp / vi) * dS)

                factor3 = (Sp + psiP) * np.linalg.inv(factor1)

                Mp = -np.linalg.inv(Sp + psiP + ((dS * (factor1)) *
                                                 np.linalg.inv(factor2)))

                Mi = -np.linalg.inv(dS + (factor3 * factor2))

                Adil1 = Id + (Sp * Mp) + (dS * Mi)

                Adil2 = Id + (Si * Mi) + (vp / vi) * (dS * (Mp - Mi))

                Ceff = (Cm)

                Ceff = self.orientational_average(
                    number_theta, number_phi, theta_max, phi_max, Ceff, 0, 1)

        else:

            number_theta = 20

            number_phi = 20

            theta_max = np.pi / 2

            phi_max = 2 * np.pi

            if vp == 0:

                Ceff = Cm

            else:

                Id = np.eye(6)

                vm = 1 - vp - vi

                psiP = (Cm * np.linalg.inv(Cp - Cm))

                psiI = (Cm * np.linalg.inv(Ci - Cm))

                dS = (Si - Sp)

                factor1 = (Sp + psiP - (vp / vi) * dS)

                factor2 = (Sp + psiI - (vp / vi) * dS)

                factor3 = ((Sp + psiP) * np.linalg.inv(factor1))

                Mp = -np.linalg.inv(Sp + psiP + ((dS * (factor1)) *
                                                 np.linalg.inv(factor2)))

                Mi = -np.linalg.inv(dS + (factor3 * factor2))

                Adil1 = Id + (Sp * Mp) + (dS * Mi)

                Adil2 = Id + (Si * Mi) + (vp / vi) * (dS * (Mp - Mi))

                Ceff = ((vm * Cm + vi * (Ci * Adil2) + vp * (Cp * Adil1)) *
                        np.linalg.inv(vm * Id + vi * Adil2 + vp * Adil1))

                Ceff = self.orientational_average(
                    number_theta, number_phi, theta_max, phi_max, Ceff, 0, 1)

        return Ceff

    def EshelbyInt(self, a, nu):

        # Compute the Eshelby tensor for general ellipsoidal inclusions.
        # Source: Mura(1987) p74-84
        #
        # a : a vector of length three containing the aspect ratios of the
        #    ellipsoidal inclusions.
        # nu : Poisson's ratio of the (isotropic) matrix material.

        A, B = self.computeEshelbyIntegrals(a)

        a = np.power(a, 2)

        c_offd = 1 / (8 * np.pi * (1 - nu))  # % Off-diagonal coefficient

        c_diag = 3 * c_offd  # % Diagonal coefficients

        c_maj = (1 - 2 * nu) * c_offd

        S = np.zeros([6, 6])

        S[0, 0] = c_maj * A[0] + c_diag * a[0] * B[0, 0]

        S[1, 1] = c_maj * A[1] + c_diag * a[1] * B[1, 1]

        S[2, 2] = c_maj * A[2] + c_diag * a[2] * B[2, 2]

        S[0, 1] = -c_maj * A[0] + c_offd * a[1] * B[0, 1]

        S[0, 2] = -c_maj * A[0] + c_offd * a[2] * B[0, 2]

        S[1, 0] = -c_maj * A[1] + c_offd * a[0] * B[1, 0]

        S[1, 2] = -c_maj * A[1] + c_offd * a[2] * B[1, 2]

        S[2, 0] = -c_maj * A[2] + c_offd * a[0] * B[2, 0]

        S[2, 1] = -c_maj * A[2] + c_offd * a[1] * B[2, 1]

        c_offd = c_offd * 0.5

        c_maj = c_maj * 0.5

        S[3, 3] = (a[1] + a[2]) * B[1, 2] * c_offd + c_maj * (A[1] + A[2])

        S[4, 4] = (a[0] + a[2]) * B[0, 2] * c_offd + c_maj * (A[0] + A[2])

        S[5, 5] = (a[0] + a[1]) * B[0, 1] * c_offd + c_maj * (A[0] + A[1])

        # Matrix notation

        S[3, :] = S[3, :] * 2

        S[4, :] = S[4, :] * 2

        S[5, :] = S[5, :] * 2

        return np.mat(S)

    def Rdsimpson13(self, f, a, b, n):
        # n debe ser par
        # valor del intervalo equiespaciado
        h = (b - a) / n
        # se discretiza  el intervalo en tamaños de h
        x = np.linspace(a, b, n + 1)
        # x=a:h:b;
        # Se crean vectores con los valores pares e inpares de x
        # para la integral
        par = x[1:-1:2]
        impar = x[2:-2:2]
        # Se calcula la integral

        Integral = (h / 3) * (f(x[0]) + 4 * np.sum(f(par)
                                                   ) + 2 * np.sum(f(impar)) + f(x[-1]))
        return Integral

    def impro_integrationv2(self, f, n):
        # Integrate fun_to_int from x=0 to x=cutoff, then from x=cutoff to
        # x=infinity.
        cutoff = 4.0

        def f2(y): return f(y**-1) * (y**-2)

        Int1 = self.Rdsimpson13(f, 0, cutoff, n)

        Int2 = self.Rdsimpson13(f2, 1e-9, 1.0 / cutoff, n)

        integral = Int1 + Int2

        return integral

    def computeEshelbyIntegrals(self, a):
        # Numerically compute the integrals I1, I2, I3, I11, I12, ... I33 required
        # for computing Eshelby's tensor. The computation is based on Mura(1987),
        # This function provides a straightforward (inefficient) implementation.
        # Its advantage is that it works for arbitrary strictly positive aspect
        # ratios.
        #
        # A: 3x1 vector containing integrals I_i.
        # B: 3x3 matrix containing integrals I_ij.
        def FD(s): return (
            np.sqrt((a[0]**2 + s) * (a[1]**2 + s) * (a[2]**2 + s)))

        def F1(s): return (1 / ((a[0]**2 + s) * FD(s)))
        def F2(s): return (1 / ((a[1]**2 + s) * FD(s)))
        # F3 = lambda s: 1/((a(2).^2+s).*FD(s))
        def F11(s): return (1 / ((a[0]**2 + s) * (a[0]**2 + s) * FD(s)))
        def F12(s): return (1 / ((a[0]**2 + s) * (a[1]**2 + s) * FD(s)))
        def F13(s): return (1 / ((a[0]**2 + s) * (a[2]**2 + s) * FD(s)))
        # F21 = lambda s: 1./((a[1]**2+s)*(a[0]**2+s)*FD(s))
        def F22(s): return (1 / ((a[1]**2 + s) * (a[1]**2 + s) * FD(s)))
        def F23(s): return (1 / ((a[1]**2 + s) * (a[2]**2 + s) * FD(s)))
        # F31 = lambda s: 1./((a(3).^2+s).*(a(1).^2+s).*FD(s))
        # F32 = lambda s: 1./((a(3).^2+s).*(a(2).^2+s).*FD(s))
        def F33(s): return (1 / ((a[2]**2 + s) * (a[2]**2 + s) * FD(s)))
        c = 2 * np.pi * np.prod(a)
        A = np.zeros([3, 1])
        B = np.zeros([3, 3])
        n = 1000000

        A[0, 0] = c * self.impro_integrationv2(F1, n)
        A[1, 0] = c * self.impro_integrationv2(F2, n)
        A[2, 0] = 4 * np.pi - A[0] - A[1]  # quadgk(F3, 0, Inf);
        # B = c*[ quadgk(F11, 0, Inf) quadgk(F12, 0, Inf) quadgk(F13, 0, Inf);
        #      quadgk(F21, 0, Inf) quadgk(F22, 0, Inf) quadgk(F23, 0, Inf);
        #      quadgk(F31, 0, Inf) quadgk(F32, 0, Inf) quadgk(F33, 0, Inf)];
        B[0, 0] = c * self.impro_integrationv2(F11, n)
        B[1, 1] = c * self.impro_integrationv2(F22, n)
        B[2, 2] = c * self.impro_integrationv2(F33, n)
        B[0, 1] = c * self.impro_integrationv2(F12, n)
        B[0, 2] = c * self.impro_integrationv2(F13, n)
        B[1, 2] = c * self.impro_integrationv2(F23, n)
        B[1, 0] = B[0, 1]
        B[2, 0] = B[0, 2]
        B[2, 1] = B[1, 2]
        return np.mat(A), np.mat(B)

    def Isotropic(self, Ey, poism):
        # Compute constitutive tensor of an isotropic material
        # INPUT
        # Ey : Young's Modulus
        # poism : Poisson's ration
        # OUTPUT
        # C: Constitutive isotropic tensor 6x6

        G = Ey / (2 * (1 + poism))

        S = np.zeros([6, 6])

        S[0, 0] = 1 / Ey

        S[1, 1] = 1 / Ey

        S[2, 2] = 1 / Ey

        S[3, 3] = 1 / G

        S[4, 4] = 1 / G

        S[5, 5] = 1 / G

        S[0, 1] = -poism / Ey

        S[0, 2] = S[0, 1]

        S[1, 2] = S[0, 1]

        S[1, 0] = S[0, 1]

        S[2, 0] = S[0, 2]

        S[2, 1] = S[1, 2]

        return np.mat(np.linalg.inv(S))

    def ComputeMechanicalProps(self):
        # Tensors
        Cm = self.Isotropic(self.EMatrix, self.NuMatrix)

        Cp = self.Isotropic(self.ECnt, self.NuCnt)

        # Interphase

        Kappa = self.Lcnt / self.Dcnt

        Deq = self.Dcnt * Kappa**(1. / 3.)

        Ci = self.Isotropic(self.EInt, self.NuInt)

        Sp = self.EshelbyInt([1, Kappa, 1], self.NuMatrix)

        Si = self.EshelbyInt([1, Kappa, 1], self.NuMatrix)

        lambd = self.Intert / Deq

        fi = self.softinterphase(Kappa, lambd, self.fc)

        Ceff = self.ellipsoidalinter_random2(Cp, Ci, Cm, Si, Sp, self.fc, fi)

        E, nu = self.computeEngineeringConstantsSqrt2(Ceff)

        return E[0], nu[0]


if __name__ == "__main__":

    # Matrix
    EMatrix = 2.5e9             # Young´s modulus MPa

    NuMatrix = 0.28               # Poisson´s ratio

    # MWCNTs

    # Isotropic material
    ECnt = 700e9                # Young´s modulus MPa

    NuCnt = 0.3                  # Poisson´s ratio

    # Mechanical properties
    Lcnt = 3.20995854347668e-6                # Length [microns]

    Dcnt = 10.3538294770233e-9              # Diameter [nm]

    EInt = 2.17e9                # Young´s modulus MPa

    NuInt = NuMatrix                  # Poisson´s ratio

    Intert = 31e-9                      # Interphase thickness nm

    type_inter = 1               # 1- Soft interphase, 2- Hard interphase

    fc = 0.01

    HMethod = MFH(Lcnt, Dcnt, EMatrix, NuMatrix, ECnt,
                  NuCnt, EInt, NuInt, Intert, fc)

    E, nu = HMethod.ComputeMechanicalProps()
    


