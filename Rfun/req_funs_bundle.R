
"albedo_mod" <-
function (trans, sunalt, albedot1) 
{
    x <- seq(0, 90, 2)
    y <- seq(0, 1, 0.05)
    xx <- sort(rep(x, length(y)))
    yy <- rep(y, length(x))
    alb <- rep(NA, length(trans))
    k <- which(sunalt > 0 & is.finite(trans) & trans <= 1.01)
    k <- k[order(sunalt[k])]
    dip <- diag(interp(xx, yy, albedot1, xo = sunalt[k], yo = trans[k])$z)
    alb[sort(k)] <- dip[order(k)]
    return(alb)
}
"blwhf" <-
function (Ts, Ta, rh, Fc) 
{
    ew <- 10^((0.7859 + 0.03477 * Ta)/(1 + 0.00412 * Ta))
    fw <- 1 + 1e-06 * P_default * (4.5 + 6e-04 * Ta^2)
    ew <- fw * ew
    rw <- eps_air * ew/(P_default - ew)
    r <- (rh/100) * rw
    e_a <- r * P_default/(eps_air + r)
    ta <- Ta + CtoK
    ts <- Ts + CtoK
    lwest <- -emiss_lw * sigmaSB * ta^4 * (0.39 - 0.05 * sqrt(e_a)) * 
        Fc - 4 * emiss_lw * sigmaSB * ta^3 * (Ts - Ta)
    return(lwest)
}
"cdntc" <-
function (sp, z, Ta) 
{
    tol <- 1e-05
    visc <- 1.326e-05 * (1 + 0.006542 * Ta + 8.301e-06 * Ta^2 - 
        4.84e-09 * Ta^3)
    sp[which(sp == 0)] <- 0.1
    ustaro <- array(0, dim(sp))
    ustarn <- 0.036 * sp
    iter <- 0
    while (abs(ustarn - ustaro) > tol) {
        ustaro <- ustarn
        z0 <- Charnock_alpha * ustaro^2/g + R_roughness * visc/ustaro
        ustarn <- sp * (kappa/log(z/z0))
        iter <- iter + 1
    }
    sqrcd <- kappa/log((10)/z0)
    cd <- sqrcd^2
    u10 = ustarn/sqrcd
    return(list(cd, u10))
}
"cdnve" <-
function (sp, z) 
{
    A <- 0.002717
    B <- 0.000142
    C <- 7.64e-05
    a <- log(z/10)/kappa
    tol <- 0.001
    u10o <- rep(0, length(sp)) + 0.1
    cd <- (A/u10o + B + C * u10o)
    u10 <- sp/(1 + a * sqrt(cd))
    ii <- abs(u10 - u10o) > tol
    while (any(ii == TRUE)) {
        u10o <- u10
        cd <- (A/u10o + B + C * u10o)
        u10 <- sp/(1 + a * sqrt(cd))
        ii <- abs(u10 - u10o) > tol
    }
    return(list(cd, u10))
}
"cloudcor" <-
function (C, lat) 
{
    ls1 <- c(0, 7, 15, 25, 35, 45, 55, 65, 75)
    ls2 <- c(0.5, 0.52, 0.59, 0.63, 0.68, 0.72, 0.76, 0.8, 0.84)
    a1 <- ls2[which.max(1/(-ls1 + lat))]
    Fc <- 1 - a1 * C
    return(Fc)
}
"convection_v12" <-
function (Tz_in, Cz_in, Sz_in, Pz_in, Chlz_in, PPz_in, DOPz_in, 
    DOCz_in, Tprof_prev, Vz, Cw, f_par, lambdaz_wtot_avg, zz, 
    swa_b0, tracer_switch, springautumn) 
{
    Trhomax <- 3.98
    Nz <- length(zz)
    dz <- zz[2] - zz[1]
    Rho <- rho(0, Tz_in, 0)
    d_rho <- c(diff(Rho), 1)
    while (length(which(d_rho < 0))) {
        blnUnstb_layers <- 1 * (d_rho <= 0)
        A_Unstb <- which(diff(c(0, blnUnstb_layers)) == 1)
        B_Unstb <- which(diff(c(0, blnUnstb_layers)) == -1) - 
            1
        for (n in 1:length(A_Unstb)) {
            j <- c(A_Unstb[n]:(B_Unstb[n] + 1))
            Tz_in[j] <- weighted.mean(Tz_in[j], Vz[j])
            if (tracer.switch == 1) 
                Cz_in[j] <- weighted.mean(Cz_in[j], Vz[j])
            Sz_in[j] <- weighted.mean(Sz_in[j], Vz[j])
            Pz_in[j] <- weighted.mean(Pz_in[j], Vz[j])
            Chlz_in[j] <- weighted.mean(Chlz_in[j], Vz[j])
            PPz_in[j] <- weighted.mean(PPz_in[j], Vz[j])
            DOPz_in[j] <- weighted.mean(DOPz_in[j], Vz[j])
            DOCz_in[j] <- weighted.mean(DOCz_in[j], Vz[j])
        }
        Rho <- rho(0, Tz_in, 0)
        d_rho <- c(diff(Rho), 1)
    }
    if (springautumn == 1) {
        jumpinx <- which(((Tprof.prev >= Trhomax) & (Tz_in < 
            Trhomax)) | ((Tprof.prev <= Trhomax) & (Tz_in > Trhomax)))
        if (length(jumpinx) > 0) {
            sellay <- jumpinx[1]:length(Tz_in)
            intSign <- sign(Trhomax - Tz_in[jumpinx[1]])
            XE_turn <- cumsum((Tz_in[sellay] - Trhomax) * Vz[sellay] * 
                Cw * intSign)
            Dummy <- which(XE_turn > 0)
            if (length(Dummy) == 0) {
                Tz_in[sellay] <- Trhomax
                if (intSign == 1) {
                  Tz_in[jumpinx[1]] <- Tz_in[jumpinx[1]] + intSign * 
                    XE_turn[length(XE_turn)]/(Vz[jumpinx[1]] * 
                    Cw)
                }
                else {
                  distrib_key <- (f.par * exp(c(0, -lambdaz.wtot.avg) * 
                    c(zz, zz[length(zz)] + dz)) + (1 - f.par) * 
                    exp(c(0, rep(-swa.b0, l = Nz)) * c(zz, max(zz) + 
                      dz)))
                  Tz_in[sellay] <- Tz_in[sellay] + (-diff(intSign * 
                    XE_turn[length(XE_turn)] * (distrib_key[jumpinx[1]:length(distrib_key)]/distrib_key[jumpinx[1]]))/(Vz[sellay] * 
                    Cw))
                }
            }
            else {
                Tz_in[jumpinx[1]:(jumpinx[1] + Dummy[1] - 1 - 
                  1)] <- Trhomax
                Tz_in[jumpinx[1] + Dummy[1] - 1] <- Trhomax + 
                  intSign * (XE_turn[Dummy[1]])/(Vz[jumpinx[1] + 
                    Dummy[1] - 1] * Cw)
            }
        }
    }
    Tz <- Tz_in
    Cz <- Cz_in
    Sz <- Sz_in
    Pz <- Pz_in
    PPz <- PPz_in
    Chlz <- Chlz_in
    DOPz <- DOPz_in
    DOCz <- DOCz_in
    return(list(Tz, Cz, Sz, Pz, Chlz, PPz, DOPz, DOCz))
}
"d2r" <-
function (deg) 
{
    (pi/180) * deg
}
"gIyear" <-
function (lat, z, daynr = 1:365, hseq = seq(0, 24, 1/6), Tlk = NULL) 
{
    if (length(lat) != 1 || length(z) != 1) 
        stop("Expected lat and z to be of length 1.")
    if (is.null(Tlk)) {
        message("--- Using the default monthly Tlk values. ---")
        message("(2.2,2.4,3.0,3.2,3.5,3.6,3.8,3.9,3.5,3.0,2.4,2.2)")
        Tlk <- c(2.2, 2.4, 3, 3.2, 3.5, 3.6, 3.8, 3.9, 3.5, 3, 
            2.4, 2.2)
        year <- as.numeric(format(Sys.Date(), "%Y"))
        x <- seq(ISOdate(year, 1, 15), ISOdate(year, 12, 15), 
            by = "month")
        x <- as.numeric(x - ISOdate(year, 1, 1))
        Tlk.spline <- splinefun(x, Tlk, method = "periodic")
        Tlk <- Tlk.spline(daynr)
        rm(year, x, Tlk.spline)
    }
    j <- (2 * pi * daynr)/365.25
    e <- 1 + 0.03344 * cos(j - 0.048869)
    sI <- 1367 * e
    declin <- asin(0.3978 * sin(j - 1.4 + (0.0355 * sin(j - 0.0489))))
    rm(e, j)
    ha <- 0.261799 * (hseq - 12)
    ha <- matrix(ha, nrow = length(daynr), ncol = length(hseq), 
        byrow = TRUE)
    sh <- cos(lat) * cos(declin) * cos(ha) + sin(lat) * sin(declin)
    ash <- asin(sh)
    dh0 <- 0.061359 * ((0.1594 + 1.123 * ash + 0.065656 * ash^2)/(1 + 
        28.9344 * ash + 277.3971 * ash^2))
    h0ref <- ash + dh0
    m <- exp(-z/8434.5)/(sin(h0ref) + 0.50572 * (h0ref + 6.07995)^-1.6364)
    hrows <- which(m > 20)
    lrows <- which(m <= 20)
    ot <- matrix(nrow = dim(m)[1], ncol = dim(m)[2])
    ot[lrows] <- 1/(6.6296 + 1.7513 * m[lrows] - 0.1202 * m[lrows]^2 + 
        0.0065 * m[lrows]^3 - 0.00013 * m[lrows]^4)
    ot[hrows] <- 1/(10.4 + 0.718 * m[hrows])
    bI <- sI * exp(-0.8662 * Tlk * m * ot) * sh
    bI[bI < 0] <- 0
    tn <- -0.015843 + 0.030543 * Tlk + 0.0003797 * Tlk^2
    A1 <- 0.26463 - 0.061581 * Tlk + 0.0031408 * Tlk^2
    x <- A1 * tn
    A1[x < 0.0022] <- 0.0022/tn[x < 0.0022]
    A2 <- 2.0402 + 0.018945 * Tlk - 0.011161 * Tlk^2
    A3 <- -1.3025 + 0.039231 * Tlk + 0.0085079 * Tlk^2
    fd <- A1 + A2 * sh + A3 * sh^2
    rm(A2, A3)
    dI <- sI * tn * fd
    dI[dI < 0] <- 0
    gI <- bI + dI
    SoAl <- 60 * asin(sh)
    return(list(sh, gI))
}
"heatflux_v11" <-
function (vec) 
{
    A <- hfbulktc_speed(vec[7], 10, vec[4], 2, vec[5], 2, vec[6], 
        Tw)
    Qsl <- A[[1]] + A[[2]]
    Qlw <- blwhf(Tw, vec[4], vec[5], cloudcor(vec[3], lat))
    dv <- chron(731276, origin = c(month = 12, day = 31, year = 1999))
    dv <- chron(vec[1], origin = c(month = 12, day = 31, year = 1999))
    yr <- dv[1]
    resol_t <- 24 * 60/5
    obj <- soradna1(dv, resol_t, 0, lat)
    alt <- obj[[1]]
    rad <- obj[[2]]
    inx <- which(alt < 0)
    alt[inx] <- 0
    Dayfrac <- 1 - length(inx)/length(alt)
    alt_trsh <- 15
    inx2 <- which(alt > alt_trsh)
    Dayfracheating <- length(inx2)/length(alt)
    if (vec[3] < 0.3) {
        CloudCorr <- 1
    }
    else {
        CloudCorr <- (1 - 0.62 * vec[3] + 0.0019 * max(alt))
    }
    if (is.na(vec[2])) {
        Trnsmiss <- rep((0.377 + 0.00513 * max(alt)) * CloudCorr, 
            length = length(alt))
    }
    else {
        Qmm <- mean(rad)
        Qmo <- (1e+06/(24 * 60 * 60)) * vec[2]
        Trnsmiss <- rep(min((Qmo/Qmm), 1), length(alt))
    }
    if (WEQs > 0) {
        alb <- alb_melt_snow
    }
    else {
        if (Hi > 0) {
            alb <- alb_melt_ice
        }
        else alb <- albedo_mod(Trnsmiss, alt, albedot1)
    }
    Qma <- (1 - alb) * rad * Trnsmiss
    Qma[which(is.na(Qma))] <- 0
    Qsw <- mean(Qma)
    if (vec[7] > 0) {
        obj <- cdnve(vec[7], 10)
        cd <- obj[[1]]
        u10 <- obj[[2]]
        tau <- rho_air * (cd * u10^2)
    }
    else {
        tau <- 0
    }
    return(c(Qsw, Qlw, Qsl, tau, Dayfrac, Dayfracheating))
}
"hfbulktc_speed" <-
function (ur, zr, Ta, zt, rh, zq, Pa, Ts) 
{
    M <- length(ur)
    tol <- 0.001
    onethird <- 1/3
    o61 <- 1/eps_air - 1
    visc <- 1.326e-05 * (1 + 0.006542 * Ta + 8.301e-06 * Ta^2 - 
        4.84e-09 * Ta^3)
    Qsats <- Qsat_coeff * qsat(Ts, Pa)
    Q <- (0.01 * rh) * qsat(Ta, Pa)
    T <- Ta + CtoK
    Tv <- T * (1 + o61 * Q)
    rho <- (100 * Pa)/(gas_const_R * Tv)
    Dt <- (Ta + 0.0098 * zt) - Ts
    Dq <- Q - Qsats
    S <- sqrt(ur^2 + min_gustiness^2)
    cdnhf <- sqrt((cdntc(S, zr, Ta))[[1]])
    z0t <- 7.5 * 10^(-5)
    ctnhf <- kappa/log(zt/z0t)
    z0q <- z0t
    cqnhf <- kappa/log(zq/z0q)
    U_star <- cdnhf * S
    T_star <- ctnhf * Dt
    Q_star <- cqnhf * Dq
    Dter <- 0
    Dqer <- 0
    Reu = 0
    Ret = 0
    Req = 0
    cvrgnce <- 1
    ii <- 0
    while (cvrgnce) {
        ReuO <- Reu
        RetO <- Ret
        ReqO <- Req
        bs <- g * (T_star * (1 + o61 * Q) + o61 * T * Q_star)/Tv
        L <- (U_star^2)/(kappa * bs)
        index_limit <- (L < zr/3 & L > 0)
        L[index_limit] = zr[index_limit]/3
        zetu <- zr/L
        zett <- zt/L
        zetq <- zq/L
        z0 <- (Charnock_alpha/g) * U_star^2 + R_roughness * visc/U_star
        cdnhf <- kappa/(log(zr/z0) - psi(zetu, "utc"))
        U_star <- cdnhf * S
        Reu <- z0 * U_star/visc
        Ret <- LKB(Reu)[[1]]
        Req <- LKB(Reu)[[2]]
        z0t <- visc * Ret/U_star
        z0q <- visc * Req/U_star
        cthf <- kappa/(log(zt/z0t) - psi(zett, "ttc"))
        cqhf <- kappa/(log(zq/z0q) - psi(zetq, "ttc"))
        T_star <- cthf * (Dt + Dter)
        Q_star <- cqhf * (Dq + Dqer)
        Ws <- U_star * (-CVB_depth/(kappa * L))^onethird
        wg <- rep(min_gustiness, length = M)
        j <- which(zetu < 0)
        wg[j] <- max(min_gustiness, beta_conv * Ws[j])
        S <- sqrt(ur^2 + wg^2)
        cvrgnce <- abs(Reu - ReuO) > tol | abs(Ret - RetO) > 
            tol | abs(Req - ReqO) > tol
        ii = ii + 1
        if (ii > 80) {
            warning("Iteration did not converge (hfbulktc_speed)!")
            break
        }
    }
    Le <- (2.501 - 0.00237 * (Ts - Dter)) * 10^6
    Hs <- rho * cp * U_star * T_star
    Hl <- rho * Le * U_star * Q_star
    CD <- (U_star/S)^2
    CT <- U_star * T_star/(S * (Dt + Dter))
    CQ <- U_star * Q_star/(S * (Dq + Dqer))
    stress <- rho * CD * S * ur
    RI <- g * zr * ((Dt + Dter) + o61 * T * (Dq + Dqer))/(Tv * 
        S^2)
    W <- 1.61 * U_star * Q_star + (1 + 1.61 * Q) * U_star * T_star/T
    Hl_webb <- rho * Le * W * Q
    A <- list(Hs, Hl, Hl_webb, stress, U_star, T_star, Q_star, 
        L, zetu, CD, CT, CQ, RI)
    return(A)
}
"LKB" <-
function (Reu) 
{
    Ret <- Reu
    Ret[] <- 1
    Req <- 0.292 * Ret[]
    j <- which(Reu > 0.11 & Reu <= 0.825)
    Ret[j] <- 1.376 * Reu[j]^0.929
    Req[j] <- 1.808 * Reu[j]^0.826
    j <- which(Reu > 0.825 & Reu <= 3)
    Ret[j] <- 1.026/Reu[j]^0.599
    Req[j] <- 1.393/Reu[j]^0.528
    j <- which(Reu > 3 & Reu <= 10)
    Ret[j] <- 1.625/Reu[j]^1.018
    Req[j] <- 1.956/Reu[j]^0.87
    j <- which(Reu > 10 & Reu <= 30)
    Ret[j] <- 4.661/Reu[j]^1.475
    Req[j] <- 4.994/Reu[j]^1.297
    j <- which(Reu > 30)
    Ret[j] <- 34.904/Reu[j]^2.067
    Req[j] <- 30.79/Reu[j]^1.845
    return(list(Ret, Req))
}
"Ppart" <-
function (vf, TIP, Psat, Fmax, rho_sed, Fstable) 
{
    N <- length(TIP)
    Pdiss <- rep(NA, l = N)
    for (w in 1:N) {
        a <- vf[w] - 1
        b <- TIP[w] + (vf[w] - 1) * Psat - vf[w] * rho_sed * 
            (Fmax + Fstable)
        c <- Psat * TIP[w] - vf[w] * rho_sed * Fstable * Psat
        Pdiss[w] <- max(as.numeric(1/polyroot(c(a, b, c))))
    }
    Pfpart = (TIP - (1 - vf) * Pdiss)/(rho_sed * vf)
    return(list(Pdiss, Pfpart))
}
"psi" <-
function (zet, meth) 
{
    c13 <- 1/3
    sq3 <- sqrt(3)
    y <- -4.7 * zet
    j <- which(zet < 0)
    zneg <- zet[j]
    x <- (1 - 16 * zneg)^0.25
    if (meth == "ttc") 
        y1 <- 2 * log((1 + x^2)/2)
    if (meth == "utc") 
        y1 <- 2 * log((1 + x)/2) + log((1 + x^2)/2) - 2 * atan(x) + 
            pi/2
    x <- (1 - 12.87 * zneg)^c13
    y2 <- 1.5 * log((x^2 + x + 1)/3) - sq3 * atan((2 * x + 1)/sq3) + 
        pi/sq3
    F <- 1/(1 + zneg^2)
    y[j] <- F * y1 + (1 - F) * y2
    return(y)
}
"qsat" <-
function (Ta, Pa) 
{
    ew = 6.1121 * (1.0007 + 3.46e-06 * Pa) * exp((17.502 * Ta)/(240.97 + 
        Ta))
    q = 0.62197 * (ew/(Pa - 0.378 * ew))
    return(q)
}
"r2d" <-
function (rad) 
{
    (180/pi) * rad
}
"rho" <-
function (S = 35, T = 25, P = 0) 
{
    rhow = 999.842594 + 0.06793952 * T - 0.00909529 * T^2 + 0.0001001685 * 
        T^3 - 1.120083e-06 * T^4 + 6.536332e-09 * T^5
    A = 0.824493 - 0.0040899 * T + 7.6438e-05 * T^2 - 8.2467e-07 * 
        T^3 + 5.3875e-09 * T^4
    B = -0.00572466 + 0.00010227 * T - 1.6546e-06 * T^2
    C = 0.00048314
    rho0 = rhow + A * S + B * S^(3/2) + C * S^2
    Ksbmw = 19652.21 + 148.4206 * T - 2.327105 * T^2 + 0.01360477 * 
        T^3 - 5.155288e-05 * T^4
    Ksbm0 = Ksbmw + S * (54.6746 - 0.603459 * T + 0.0109987 * 
        T^2 - 6.167e-05 * T^3) + S^(3/2) * (0.07944 + 0.016483 * 
        T - 0.00053009 * T^2)
    Ksbm = Ksbm0 + P * (3.239908 + 0.00143713 * T + 0.000116092 * 
        T^2 - 5.77905e-07 * T^3) + P * S * (0.0022838 - 1.0981e-05 * 
        T - 1.6078e-06 * T^2) + P * S^(3/2) * 0.000191075 + P * 
        P * (8.50935e-05 - 6.12293e-06 * T + 5.2787e-08 * T^2) + 
        P^2 * S * (-9.9348e-07 + 2.0816e-08 * T + 9.1697e-10 * 
            T^2)
    rho = rho0/(1 - P/Ksbm)
    return(rho)
}
"soradna1" <-
function (dat, res, long, lat) 
{
    dat <- dat + seq(0, 1, 1/res)[-(res + 1)]
    SC <- seconds(dat)
    MN <- minutes(dat)
    HR <- hours(dat)
    D <- M <- Y <- rep(NA, length = length(MN))
    D[] <- as.numeric(unlist(strsplit(as.character(dates(dat)[1]), 
        "/"))[2])
    M[] <- as.numeric(unlist(strsplit(as.character(dates(dat)[1]), 
        "/"))[1])
    Y[] <- as.numeric(unlist(strsplit(as.character(dates(dat)[1]), 
        "/"))[3])
    if (Y[1] < 1900) 
        Y[] <- Y[1] + 2000
    LONG <- rep(long, length = length(SC))
    LAT <- rep(lat, length = length(SC))
    DTR <- 3.14159265/180
    RTD <- 1/DTR
    UT <- HR + (MN + SC/60)/60
    JED <- 367 * Y - floor(7 * (Y + floor((M + 9)/12))/4) + floor(275 * 
        M/9) + D + 1721013 + UT/24
    T <- (JED - 2415020)/36525
    G <- 358.475833 + 35999.04975 * T - 0.00015 * T^2
    NG <- floor(G/360)
    G <- (G - NG * 360) * DTR
    L <- 279.696678 + 36000.76892 * T + 0.000303 * T^2
    NL <- floor(L/360)
    L <- (L - NL * 360) * DTR
    JUP <- 225.444651 + 2880 * T + 154.906654 * T
    NJUP <- floor(JUP/360)
    JUP <- (JUP - NJUP * 360) * DTR
    NM <- 259.183275 - 1800 * T - 134.142008 * T + 0.002078 * 
        T^2
    NNM <- floor(NM/360)
    NM <- (NM - NNM * 360 + 360) * DTR
    V <- 212.603219 + 58320 * T + 197.803875 * T + 0.001286 * 
        T^2
    NV <- floor(V/360)
    V <- (V - NV * 360) * DTR
    THETA <- 0.39793 * sin(L) + 0.009999 * sin(G - L) + 0.003334 * 
        sin(G + L) - 0.000208 * T * sin(L) + 4.2e-05 * sin(2 * 
        G + L) - 4e-05 * cos(L) - 3.9e-05 * sin(NM - L) - 3e-05 * 
        T * sin(G - L) - 1.4e-05 * sin(2 * G - L) - 1e-05 * cos(G - 
        L - JUP) - 1e-05 * T * sin(G + L)
    RHO <- 1.000421 - 0.033503 * cos(G) - 0.00014 * cos(2 * G) + 
        8.4e-05 * T * cos(G) - 3.3e-05 * sin(G - JUP) + 2.7e-05 * 
        sin(2 * G - 2 * V)
    DECL <- asin(THETA/sqrt(RHO))
    L <- 276.697 + 0.98564734 * (JED - 2415020)
    L <- (L - 360 * floor(L/360)) * DTR
    EQT <- -97.8 * sin(L) - 431.3 * cos(L) + 596.6 * sin(2 * 
        L) - 1.9 * cos(2 * L) + 4 * sin(3 * L) + 19.3 * cos(3 * 
        L) - 12.7 * sin(4 * L)
    EQT <- EQT/60
    L <- L * RTD
    GHA <- 15 * (UT - 12) + 15 * EQT/60
    LHA <- GHA - LONG
    RV <- sqrt(RHO)
    SZ <- sin(DTR * LAT) * sin(DECL) + cos(DTR * LAT) * cos(DECL) * 
        cos(DTR * LHA)
    z <- RTD * asin(SZ)
    sorad <- rep(0, length = length(z))
    ii <- which(z > 0)
    sorad[ii] <- (Solar_const/RV[ii]^2) * sin(DTR * z[ii])
    return(list(z, sorad))
}
