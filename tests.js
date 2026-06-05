// tests.js — pure Node tests for engine.js (DTA Meta-Analysis Suite).
// Every expected value is hand-derived independently of the implementation.
// Run: node tests.js   ->   "N passed, M failed" and exit(M===0?0:1).

const E = require('./engine.js');
const {
  logit, invLogit, wilsonCI, calcStudyStats, univariateRE,
  bivariateModel, calcSROC, thresholdEffectTest, normalCDF, deeksFunnelTest
} = E;

let passed = 0, failed = 0;
function approx(name, got, exp, tol) {
  tol = tol == null ? 1e-9 : tol;
  if (typeof got === 'number' && isFinite(got) && Math.abs(got - exp) <= tol) {
    passed++; console.log('PASS ' + name + '  (got ' + got + ')');
  } else {
    failed++; console.log('FAIL ' + name + '  got ' + got + ' expected ' + exp + ' (tol ' + tol + ')');
  }
}
function ok(name, cond, detail) {
  if (cond) { passed++; console.log('PASS ' + name); }
  else { failed++; console.log('FAIL ' + name + (detail ? '  ' + detail : '')); }
}

// ----------------------------------------------------------------------------
// 1. logit / invLogit round-trip and clamping
// ----------------------------------------------------------------------------
// logit(0.8) = ln(0.8/0.2) = ln(4) = 1.3862943611198906
approx('logit(0.8) = ln(4)', logit(0.8), Math.log(4), 1e-12);
// invLogit(0) = 1/(1+e^0) = 0.5
approx('invLogit(0) = 0.5', invLogit(0), 0.5, 1e-12);
// round trip
approx('invLogit(logit(0.37)) round-trip', invLogit(logit(0.37)), 0.37, 1e-12);

// CLAMPING: logit(0) and logit(1) must NOT be +-Infinity (clamped to [1e-4, 1-1e-4]).
// logit(0) clamps p->0.0001 => ln(0.0001/0.9999) = ln(0.00010001) = -9.210340...
const lo0 = logit(0), lo1 = logit(1);
ok('logit(0) is finite (clamped, not -Inf)', isFinite(lo0), 'got ' + lo0);
ok('logit(1) is finite (clamped, not +Inf)', isFinite(lo1), 'got ' + lo1);
// Hand: ln(0.0001 / 0.9999) = ln(0.0001) - ln(0.9999) = -9.210340372 - (-0.000100005) = -9.210240367
approx('logit(0) clamped value', lo0, Math.log(0.0001 / 0.9999), 1e-12);
ok('logit(0) == -logit(1) (symmetric clamp)', Math.abs(lo0 + lo1) < 1e-12, 'lo0=' + lo0 + ' lo1=' + lo1);

// ----------------------------------------------------------------------------
// 2. normalCDF: Phi(0)=0.5 (bug if it returns 0), and known tail values.
// ----------------------------------------------------------------------------
approx('normalCDF(0) = 0.5 (NOT 0)', normalCDF(0), 0.5, 1e-6);
// Phi(1.959963985) = 0.975 ; Phi(1.96) ~ 0.9750021
approx('normalCDF(1.96) ~ 0.975', normalCDF(1.96), 0.9750021048, 1e-4);
approx('normalCDF(-1.96) ~ 0.025', normalCDF(-1.96), 0.0249978952, 1e-4);
ok('normalCDF symmetric Phi(x)+Phi(-x)=1', Math.abs((normalCDF(1.3) + normalCDF(-1.3)) - 1) < 1e-6);

// ----------------------------------------------------------------------------
// 3. Hand-worked 2x2 study: TP=80, FP=10, FN=20, TN=90.
//    diseased=100, healthy=100, total=200.
//    Continuity-corrected (+0.5): tpC=80.5, fpC=10.5, fnC=20.5, tnC=90.5
//                                 diseasedC=101, healthyC=101
//    sens (raw)   = 80/100 = 0.8
//    spec (raw)   = 90/100 = 0.9
//    sensLogit    = logit(80.5/101) = logit(0.79702970...) = ln(0.79702970/0.20297030)
//                 = ln(3.9268292...) = 1.3678322982801554
//    specLogit    = logit(90.5/101) = logit(0.89603960...) = ln(0.89603960/0.10396040)
//                 = ln(8.6190476...) = 2.153974593542403
//    DOR (cross-product, +0.5) = (tpC*tnC)/(fpC*fnC) = (80.5*90.5)/(10.5*20.5)
//                 = 7285.25 / 215.25 = 33.84552845528455
//    logDOR       = ln(33.84552845528455) = 3.521806891822558
//    varSensLogit = 1/tpC + 1/fnC = 1/80.5 + 1/20.5 = 0.012422360 + 0.048780488 = 0.061202848
//    varSpecLogit = 1/tnC + 1/fpC = 1/90.5 + 1/10.5 = 0.011049724 + 0.095238095 = 0.106287819
//    varLogDOR    = 1/80.5+1/10.5+1/20.5+1/90.5 = 0.167490667
// ----------------------------------------------------------------------------
const s = calcStudyStats({ tp: 80, fp: 10, fn: 20, tn: 90 });
approx('study sens (raw)', s.sens, 0.8, 1e-12);
approx('study spec (raw)', s.spec, 0.9, 1e-12);
approx('study sensLogit', s.sensLogit, 1.3678322982801554, 1e-9);
approx('study specLogit', s.specLogit, 2.153974593542403, 1e-9);
approx('study DOR = (tpC*tnC)/(fpC*fnC)', s.dor, 33.84552845528455, 1e-7);
approx('study logDOR = ln(DOR)', s.logDOR, 3.521806891822558, 1e-9);
approx('study varSensLogit = 1/tpC+1/fnC', s.varSensLogit, 1 / 80.5 + 1 / 20.5, 1e-12);
approx('study varSpecLogit = 1/tnC+1/fpC', s.varSpecLogit, 1 / 90.5 + 1 / 10.5, 1e-12);
approx('study varLogDOR', s.varLogDOR, 1 / 80.5 + 1 / 10.5 + 1 / 20.5 + 1 / 90.5, 1e-12);

// ----------------------------------------------------------------------------
// 4. DOR SIGN CHECK (the classic bug):  DOR = exp(mu1 + mu2), NOT exp(mu1 - mu2).
//    For a single study, exp(sensLogit + specLogit) must equal the cross-product DOR,
//    and exp(sensLogit - specLogit) must NOT.
//    exp(1.3678322982801554 + 2.153974593542403) = exp(3.5218068918...) = 33.84552845528...
//    exp(1.3678322982801554 - 2.153974593542403) = exp(-0.7861423...)  = 0.45554... (wrong)
// ----------------------------------------------------------------------------
const sumLogit = s.sensLogit + s.specLogit;
const diffLogit = s.sensLogit - s.specLogit;
approx('DOR == exp(mu1+mu2)  [correct sign]', Math.exp(sumLogit), s.dor, 1e-7);
ok('DOR != exp(mu1-mu2)  [the buggy subtraction]',
   Math.abs(Math.exp(diffLogit) - s.dor) > 1.0,
   'exp(mu1-mu2)=' + Math.exp(diffLogit) + ' vs DOR=' + s.dor);

// Bivariate pooled DOR also uses the additive identity (DOR = PLR/NLR = exp(mu1+mu2)).
// Build a tiny symmetric 3-study set so we can confirm DOR = PLR/NLR algebraically.
const triStats = [
  calcStudyStats({ tp: 80, fp: 10, fn: 20, tn: 90 }),
  calcStudyStats({ tp: 85, fp: 12, fn: 15, tn: 88 }),
  calcStudyStats({ tp: 78, fp: 9,  fn: 22, tn: 91 })
];
const bm = bivariateModel(triStats);
ok('bivariateModel returns object for k=3', bm && bm.DOR && bm.PLR && bm.NLR);
approx('pooled DOR == PLR/NLR (== exp(mu1+mu2))',
       bm.DOR.estimate, bm.PLR.estimate / bm.NLR.estimate, 1e-9);
// Independent identity: PLR = sens/(1-spec), NLR = (1-sens)/spec, DOR = sens*spec/((1-spec)*(1-sens))
const ps = bm.sensitivity.estimate, pc = bm.specificity.estimate;
approx('pooled PLR = sens/(1-spec)', bm.PLR.estimate, ps / (1 - pc), 1e-9);
approx('pooled NLR = (1-sens)/spec', bm.NLR.estimate, (1 - ps) / pc, 1e-9);
approx('pooled DOR = sens*spec/((1-spec)(1-sens))',
       bm.DOR.estimate, (ps * pc) / ((1 - pc) * (1 - ps)), 1e-9);
ok('pooled DOR > 1 for a good test', bm.DOR.estimate > 1, 'DOR=' + bm.DOR.estimate);

// ----------------------------------------------------------------------------
// 5. SROC SIGN CHECK: logit(Spec) = -logit(FPR).  Spec = 1 - FPR.
//    calcSROC builds spec-logit from logit(1-fpr); confirm the sign convention.
//    For fpr=0.2: logit(1-0.2)=logit(0.8)=ln(4)=1.386294; -logit(0.2)=-ln(0.25)=1.386294. Equal.
// ----------------------------------------------------------------------------
const fpr = 0.2;
approx('logit(1-fpr) == -logit(fpr)  [SROC sign flip]',
       logit(1 - fpr), -logit(fpr), 1e-12);
const sr = calcSROC(triStats);
ok('calcSROC returns curve', sr && Array.isArray(sr.points) && sr.points.length > 0);
ok('SROC AUC in [0.5,1]', sr.auc >= 0.5 && sr.auc <= 1, 'auc=' + sr.auc);
// qStar = invLogit(alpha/2) is in (0,1)
ok('SROC qStar in (0,1)', sr.qStar > 0 && sr.qStar < 1, 'qStar=' + sr.qStar);

// ----------------------------------------------------------------------------
// 6. Bivariate PREDICTION-REGION constant: sqrt(chi2_{0.05,2}) = sqrt(5.991465) = 2.4477468.
//    This is the 2-df radius used for a 95% bivariate region (NOT the univariate z=1.96).
//    Documented constant check so any future region code uses the right multiplier.
// ----------------------------------------------------------------------------
const CHI2_2_95 = 5.991464547107979; // qchisq(0.95, df=2)
approx('sqrt(chi2_{0.05,2}) = 2.4477468 (bivariate region radius)',
       Math.sqrt(CHI2_2_95), 2.4477468307, 1e-6);
ok('bivariate radius (2.4477) != univariate z (1.96)',
   Math.abs(Math.sqrt(CHI2_2_95) - 1.96) > 0.4);

// ----------------------------------------------------------------------------
// 7. univariateRE (DerSimonian-Laird) sanity: identical values => tau2=0, estimate=value.
//    values=[1,1,1,1], variances=[0.1,0.1,0.1,0.1]:
//      fixedEst = 1 ; Q = sum w*(v-fixed)^2 = 0 ; tau2 = max(0,(0-3)/C) = 0 ; reEst = 1.
// ----------------------------------------------------------------------------
const re = univariateRE([1, 1, 1, 1], [0.1, 0.1, 0.1, 0.1]);
approx('univariateRE identical-values estimate = 1', re.estimate, 1, 1e-12);
approx('univariateRE identical-values tau2 = 0', re.tau2, 0, 1e-12);
approx('univariateRE identical-values Q = 0', re.Q, 0, 1e-12);
approx('univariateRE identical-values I2 = 0', re.I2, 0, 1e-12);
// Two-study inverse-variance weighted mean: values=[2,4], variances=[1,1] -> fixed mean = 3.
const re2 = univariateRE([2, 4], [1, 1]);
approx('univariateRE [2,4]/[1,1] estimate = 3', re2.estimate, 3, 1e-12);

// ----------------------------------------------------------------------------
// 8. wilsonCI: x=n (all successes) upper clamps to 1; symmetric center for x=5,n=10.
//    For p=0.5,n=10,z=1.96: center weight => point estimate stays 0.5 by symmetry.
// ----------------------------------------------------------------------------
const w = wilsonCI(5, 10);
ok('wilsonCI(5,10) bounds within [0,1]', w.lower >= 0 && w.upper <= 1 && w.lower < 0.5 && w.upper > 0.5,
   JSON.stringify(w));
const wAll = wilsonCI(10, 10);
approx('wilsonCI(10,10) upper = 1', wAll.upper, 1, 1e-9);
const wZero = wilsonCI(0, 10);
approx('wilsonCI(0,10) lower = 0', wZero.lower, 0, 1e-9);
const wEmpty = wilsonCI(0, 0);
ok('wilsonCI n=0 returns {0,1} (no div-by-zero)', wEmpty.lower === 0 && wEmpty.upper === 1);

// ----------------------------------------------------------------------------
// 9. thresholdEffectTest: perfectly concordant ranks => Spearman rho = +1.
//    Build studies where higher sens pairs with higher spec monotonically.
//    sens and spec arrays are co-monotone => sum d^2 = 0 => rho = 1.
// ----------------------------------------------------------------------------
const concordant = [
  calcStudyStats({ tp: 60, fp: 40, fn: 40, tn: 60 }), // sens .6 spec .6
  calcStudyStats({ tp: 70, fp: 30, fn: 30, tn: 70 }), // sens .7 spec .7
  calcStudyStats({ tp: 80, fp: 20, fn: 20, tn: 80 }), // sens .8 spec .8
  calcStudyStats({ tp: 90, fp: 10, fn: 10, tn: 90 })  // sens .9 spec .9
];
const th = thresholdEffectTest(concordant);
approx('thresholdEffectTest concordant Spearman rho = +1', th.rho, 1, 1e-12);

// ----------------------------------------------------------------------------
// 10. deeksFunnelTest: returns nulls for k<5 (guard), object for k>=5.
// ----------------------------------------------------------------------------
const dk4 = deeksFunnelTest(triStats.concat([triStats[0]])); // k=4
ok("deeksFunnelTest k<5 returns null slope", dk4.slope === null);
// Five studies with VARYING sample sizes (so 1/sqrt(ESS) has nonzero spread => finite slope).
const five = [
  calcStudyStats({ tp: 30, fp: 20, fn: 20, tn: 30 }),   // D=50  H=50  total=100  ESS=4*50*50/100=100
  calcStudyStats({ tp: 60, fp: 40, fn: 40, tn: 60 }),   // D=100 H=100 total=200  ESS=200
  calcStudyStats({ tp: 90, fp: 60, fn: 60, tn: 90 }),   // D=150 H=150 total=300  ESS=300
  calcStudyStats({ tp: 120, fp: 80, fn: 80, tn: 120 }), // D=200 H=200 total=400  ESS=400
  calcStudyStats({ tp: 150, fp: 100, fn: 100, tn: 150 })// D=250 H=250 total=500  ESS=500
];
const dk5 = deeksFunnelTest(five);
ok('deeksFunnelTest k=5 returns finite slope', isFinite(dk5.slope), 'slope=' + dk5.slope);
ok('deeksFunnelTest ESS computed for each study', dk5.ess.length === 5);
// ESS = 4*diseased*healthy/total. First study: 4*50*50/100 = 100 ; last: 4*250*250/500 = 500.
approx('Deeks ESS = 4*D*H/total (first study = 100)', dk5.ess[0], 100, 1e-9);
approx('Deeks ESS = 4*D*H/total (last study = 500)', dk5.ess[4], 500, 1e-9);

// ----------------------------------------------------------------------------
console.log('\n' + passed + ' passed, ' + failed + ' failed');
process.exit(failed === 0 ? 0 : 1);
