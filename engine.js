// engine.js — pure diagnostic-test-accuracy meta-analysis functions, extracted verbatim from DTA (index.html).
// No DOM / no Plotly. Browser uses these globals; Node consumes via module.exports.

// Logit and inverse logit
function logit(p) {
    p = Math.max(0.0001, Math.min(0.9999, p));
    return Math.log(p / (1 - p));
}

function invLogit(x) {
    return 1 / (1 + Math.exp(-x));
}

// Wilson score interval for single proportion
function wilsonCI(x, n, z = 1.96) {
    if (n === 0) return { lower: 0, upper: 1 };
    const p = x / n;
    const denom = 1 + z * z / n;
    const center = p + z * z / (2 * n);
    const margin = z * Math.sqrt((p * (1 - p) + z * z / (4 * n)) / n);
    return {
        lower: Math.max(0, (center - margin) / denom),
        upper: Math.min(1, (center + margin) / denom)
    };
}

// Calculate study-level statistics
function calcStudyStats(study) {
    const tp = study.tp, fp = study.fp, fn = study.fn, tn = study.tn;
    const diseased = tp + fn;
    const healthy = fp + tn;
    const total = diseased + healthy;
    
    const sens = diseased > 0 ? tp / diseased : 0;
    const spec = healthy > 0 ? tn / healthy : 0;
    
    // Add 0.5 continuity correction for zero cells
    const tpC = tp + 0.5, fpC = fp + 0.5, fnC = fn + 0.5, tnC = tn + 0.5;
    const diseasedC = tpC + fnC, healthyC = fpC + tnC;
    
    const sensLogit = logit(tpC / diseasedC);
    const specLogit = logit(tnC / healthyC);
    
    // Variance of logit-transformed proportions
    const varSensLogit = 1 / tpC + 1 / fnC;
    const varSpecLogit = 1 / tnC + 1 / fpC;
    
    // DOR with continuity correction
    const dor = (tpC * tnC) / (fpC * fnC);
    const logDOR = Math.log(dor);
    const varLogDOR = 1/tpC + 1/fpC + 1/fnC + 1/tnC;
    
    // Wilson CIs
    const sensCI = wilsonCI(tp, diseased);
    const specCI = wilsonCI(tn, healthy);
    
    return {
        sens, spec, sensCI, specCI,
        sensLogit, specLogit, varSensLogit, varSpecLogit,
        dor, logDOR, varLogDOR,
        diseased, healthy, total,
        tp, fp, fn, tn
    };
}

// DerSimonian-Laird random effects for univariate meta-analysis
function univariateRE(values, variances) {
    const k = values.length;
    if (k === 0) return null;
    
    // Fixed effect weights and estimate
    const weights = variances.map(v => 1 / v);
    const sumW = weights.reduce((a, b) => a + b, 0);
    const fixedEst = weights.reduce((sum, w, i) => sum + w * values[i], 0) / sumW;
    
    // Q statistic
    const Q = weights.reduce((sum, w, i) => sum + w * Math.pow(values[i] - fixedEst, 2), 0);
    
    // DerSimonian-Laird tau-squared
    const sumW2 = weights.reduce((a, b) => a + b * b, 0);
    const C = sumW - sumW2 / sumW;
    let tau2 = Math.max(0, (Q - (k - 1)) / C);
    
    // Random effects weights and estimate
    const reWeights = variances.map(v => 1 / (v + tau2));
    const sumREW = reWeights.reduce((a, b) => a + b, 0);
    const reEst = reWeights.reduce((sum, w, i) => sum + w * values[i], 0) / sumREW;
    const reVar = 1 / sumREW;
    const reSE = Math.sqrt(reVar);
    
    // I-squared
    const I2 = k > 1 ? Math.max(0, (Q - (k - 1)) / Q * 100) : 0;
    
    return {
        estimate: reEst,
        se: reSE,
        ci: [reEst - 1.96 * reSE, reEst + 1.96 * reSE],
        tau2,
        Q,
        I2,
        k
    };
}

// Bivariate random effects model (simplified REML-like estimation)
function bivariateModel(studyStats) {
    const k = studyStats.length;
    if (k < 3) return null;
    
    // Extract logit-transformed values
    const sensLogits = studyStats.map(s => s.sensLogit);
    const specLogits = studyStats.map(s => s.specLogit);
    const varSens = studyStats.map(s => s.varSensLogit);
    const varSpec = studyStats.map(s => s.varSpecLogit);
    
    // Univariate estimates first
    const sensRE = univariateRE(sensLogits, varSens);
    const specRE = univariateRE(specLogits, varSpec);
    
    // Correlation between logit(sens) and logit(spec) across studies
    const meanSens = sensLogits.reduce((a, b) => a + b, 0) / k;
    const meanSpec = specLogits.reduce((a, b) => a + b, 0) / k;
    
    let sumProd = 0, sumSensSq = 0, sumSpecSq = 0;
    for (let i = 0; i < k; i++) {
        const ds = sensLogits[i] - meanSens;
        const dc = specLogits[i] - meanSpec;
        sumProd += ds * dc;
        sumSensSq += ds * ds;
        sumSpecSq += dc * dc;
    }
    
    const correlation = sumProd / Math.sqrt(sumSensSq * sumSpecSq);
    
    // Transform back to probability scale
    const pooledSens = invLogit(sensRE.estimate);
    const pooledSpec = invLogit(specRE.estimate);
    
    // CI on probability scale
    const sensLower = invLogit(sensRE.ci[0]);
    const sensUpper = invLogit(sensRE.ci[1]);
    const specLower = invLogit(specRE.ci[0]);
    const specUpper = invLogit(specRE.ci[1]);
    
    // Likelihood ratios
    const PLR = pooledSens / (1 - pooledSpec);
    const NLR = (1 - pooledSens) / pooledSpec;
    
    // DOR
    const DOR = PLR / NLR;
    
    // Approximate CIs for LRs using delta method
    const plrSE = PLR * Math.sqrt(Math.pow(sensRE.se, 2) + Math.pow(specRE.se, 2));
    const nlrSE = NLR * Math.sqrt(Math.pow(sensRE.se, 2) + Math.pow(specRE.se, 2));
    
    // Log-scale CIs for LRs
    const logPLR = Math.log(PLR);
    const logNLR = Math.log(NLR);
    const logPLRSE = Math.sqrt(Math.pow(sensRE.se / pooledSens, 2) + Math.pow(specRE.se / (1 - pooledSpec), 2));
    const logNLRSE = Math.sqrt(Math.pow(sensRE.se / (1 - pooledSens), 2) + Math.pow(specRE.se / pooledSpec, 2));
    
    return {
        sensitivity: {
            estimate: pooledSens,
            ci: [sensLower, sensUpper],
            se: sensRE.se,
            logit: sensRE.estimate,
            tau2: sensRE.tau2,
            I2: sensRE.I2
        },
        specificity: {
            estimate: pooledSpec,
            ci: [specLower, specUpper],
            se: specRE.se,
            logit: specRE.estimate,
            tau2: specRE.tau2,
            I2: specRE.I2
        },
        PLR: {
            estimate: PLR,
            ci: [Math.exp(logPLR - 1.96 * logPLRSE), Math.exp(logPLR + 1.96 * logPLRSE)]
        },
        NLR: {
            estimate: NLR,
            ci: [Math.exp(logNLR - 1.96 * logNLRSE), Math.exp(logNLR + 1.96 * logNLRSE)]
        },
        DOR: {
            estimate: DOR,
            ci: [DOR * 0.5, DOR * 2] // Rough approximation
        },
        correlation,
        k
    };
}

// SROC curve calculation (Moses-Littenberg approach)
function calcSROC(studyStats) {
    const k = studyStats.length;
    if (k < 2) return null;
    
    // D = logit(sens) + logit(spec), S = logit(sens) - logit(spec)
    const D = studyStats.map(s => s.sensLogit + s.specLogit);
    const S = studyStats.map(s => s.sensLogit - s.specLogit);
    
    // Regression: D on S (or S on D)
    const meanD = D.reduce((a, b) => a + b, 0) / k;
    const meanS = S.reduce((a, b) => a + b, 0) / k;
    
    let sumSD = 0, sumSS = 0;
    for (let i = 0; i < k; i++) {
        sumSD += (S[i] - meanS) * (D[i] - meanD);
        sumSS += (S[i] - meanS) * (S[i] - meanS);
    }
    
    const beta = sumSS > 0 ? sumSD / sumSS : 0;
    const alpha = meanD - beta * meanS;
    
    // Generate SROC curve points
    const srocPoints = [];
    for (let fpr = 0.01; fpr <= 0.99; fpr += 0.01) {
        const specLogit = logit(1 - fpr);
        // From the symmetric model: sensLogit = (alpha + (1-beta)*specLogit) / (1+beta)
        // Simplified: when beta ≈ 0 (symmetric), sensLogit ≈ (alpha/2) + specLogit
        const sensLogit = (alpha - (1 - beta) * specLogit) / (1 + beta);
        const sens = invLogit(sensLogit);
        if (sens > 0 && sens < 1) {
            srocPoints.push({ fpr, sens });
        }
    }
    
    // AUC approximation using trapezoidal rule
    let auc = 0;
    for (let i = 1; i < srocPoints.length; i++) {
        const dx = srocPoints[i].fpr - srocPoints[i-1].fpr;
        const avgSens = (srocPoints[i].sens + srocPoints[i-1].sens) / 2;
        auc += dx * avgSens;
    }
    
    // Q* index (where sens = spec on SROC)
    const qStar = invLogit(alpha / 2);
    
    // D-prime (discriminability index)
    const dPrime = alpha / Math.sqrt(2);
    
    return {
        points: srocPoints,
        alpha,
        beta,
        auc: Math.max(0.5, Math.min(1, auc)),
        qStar,
        dPrime
    };
}

// Threshold effect test (Spearman correlation)
function thresholdEffectTest(studyStats) {
    const k = studyStats.length;
    if (k < 4) return { rho: null, p: null };
    
    const sens = studyStats.map(s => s.sens);
    const spec = studyStats.map(s => s.spec);
    
    // Rank the values
    function rank(arr) {
        const sorted = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
        const ranks = new Array(arr.length);
        sorted.forEach((item, r) => { ranks[item.i] = r + 1; });
        return ranks;
    }
    
    const rankSens = rank(sens);
    const rankSpec = rank(spec);
    
    // Spearman correlation
    let sumD2 = 0;
    for (let i = 0; i < k; i++) {
        const d = rankSens[i] - rankSpec[i];
        sumD2 += d * d;
    }
    const rho = 1 - (6 * sumD2) / (k * (k * k - 1));
    
    // Approximate p-value using t-distribution
    const t = rho * Math.sqrt((k - 2) / (1 - rho * rho));
    // Simple approximation for p-value
    const p = 2 * (1 - normalCDF(Math.abs(t) * Math.sqrt(k / 10)));
    
    return { rho, p };
}

// Normal CDF approximation
function normalCDF(x) {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x) / Math.sqrt(2);
    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return 0.5 * (1 + sign * y);
}

// Deeks' funnel plot asymmetry test
function deeksFunnelTest(studyStats) {
    const k = studyStats.length;
    if (k < 5) return { slope: null, p: null };
    
    // ESS = 4 * (TP+FN) * (FP+TN) / (TP+FP+FN+TN)
    const ess = studyStats.map(s => 4 * s.diseased * s.healthy / s.total);
    const invSqrtESS = ess.map(e => 1 / Math.sqrt(e));
    const logDOR = studyStats.map(s => s.logDOR);
    
    // Linear regression: logDOR ~ 1/sqrt(ESS)
    const meanX = invSqrtESS.reduce((a, b) => a + b, 0) / k;
    const meanY = logDOR.reduce((a, b) => a + b, 0) / k;
    
    let sumXY = 0, sumXX = 0, sumYY = 0;
    for (let i = 0; i < k; i++) {
        const dx = invSqrtESS[i] - meanX;
        const dy = logDOR[i] - meanY;
        sumXY += dx * dy;
        sumXX += dx * dx;
        sumYY += dy * dy;
    }
    
    const slope = sumXY / sumXX;
    const intercept = meanY - slope * meanX;
    
    // R-squared and standard error of slope
    const ssRes = sumYY - slope * sumXY;
    const mse = ssRes / (k - 2);
    const seSlope = Math.sqrt(mse / sumXX);
    
    // t-test for slope
    const t = slope / seSlope;
    const p = 2 * (1 - normalCDF(Math.abs(t)));
    
    return {
        slope,
        intercept,
        p,
        ess,
        invSqrtESS,
        logDOR
    };
}

if (typeof module!=='undefined'&&module.exports){ module.exports = { logit, invLogit, wilsonCI, calcStudyStats, univariateRE, bivariateModel, calcSROC, thresholdEffectTest, normalCDF, deeksFunnelTest }; }
