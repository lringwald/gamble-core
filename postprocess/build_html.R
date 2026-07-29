#!/usr/bin/env Rscript
suppressMessages({library(data.table)})
CACHE <- "output/report_work/cache"
MAPDIR <- "output/report_work/maps"
OUTHTML <- "output/report/prior_model_fit_report.html"
dir.create(dirname(OUTHTML), showWarnings=FALSE, recursive=TRUE)
b64 <- function(p){ if(!file.exists(p)) return(""); r<-readBin(p,"raw",file.info(p)$size)
  paste0("data:image/png;base64,", base64enc::base64encode(r)) }
esc <- function(s) { s<-gsub("&","&amp;",s); s<-gsub("<","&lt;",s); gsub(">","&gt;",s) }
fnum <- function(x,d=0) formatC(x, format="f", digits=d, big.mark=",")
km2_to_Mha <- function(x) x/10000
BRs <- c("AGMIP","GLOBIOM","BIOCLIMA")
D <- lapply(setNames(BRs,BRs), function(b){ f<-file.path(CACHE,paste0(b,".rds")); if(file.exists(f)) readRDS(f) else NULL })
D <- D[!sapply(D,is.null)]

tile <- function(label,val,sub="") sprintf('<div class="tile"><div class="tval">%s</div><div class="tlabel">%s</div>%s</div>',
  val, esc(label), if(nzchar(sub)) sprintf('<div class="tsub">%s</div>',sub) else "")

branch_section <- function(b){
  R<-D[[b]]; cats<-R$cats; J<-R$J
  fitm<-as.data.table(R$fitm)[order(-obs_share)]
  ct<-as.data.table(R$class_tot); ct[,pdiff:=100*(est_median-actual)/pmax(actual,1e-9)]
  ct<-ct[order(-actual)]
  conv<-R$conv
  # stat tiles
  tiles<-paste0(
    tile("grid pixels (10 km)", fnum(R$n)),
    tile("land-use classes", J),
    tile(paste0("RE regions (",esc(R$cfg$re),")"), length(R$glev)),
    tile("McFadden R² (spatial)", sprintf("%.3f", R$McF_recon), "in-sample, exact slopes"),
    tile("R-hat < 1.01 (slopes)", sprintf("%.0f%%", 100*conv$frac_rhat_lt_1.01), sprintf("%.0f%% < 1.1", 100*conv$frac_rhat_lt_1.1)),
    tile("median ESS (slopes)", fnum(conv$median_ess), sprintf("%d posterior draws", R$ndraw_used))
  )
  # calibration table
  crows<-paste(apply(fitm,1,function(r){
    oc<-as.numeric(r[["obs_share"]]); pc<-as.numeric(r[["pred_share"]]); pr<-as.numeric(r[["pixel_corr"]])
    barw<-if(is.na(pr)) 0 else max(0,pr)*100
    sprintf('<tr><td class="cl">%s</td><td class="n">%.2f%%</td><td class="n">%.2f%%</td><td class="n">%s</td><td class="barcell"><span class="bar" style="width:%.0f%%"></span></td></tr>',
      esc(r[["class"]]), 100*oc, 100*pc, if(is.na(pr)) "–" else sprintf("%.2f",pr), barw)
  }),collapse="")
  cal<-sprintf('<table class="dt"><thead><tr><th>class</th><th class="n">observed</th><th class="n">predicted</th><th class="n">pixel r</th><th>spatial skill</th></tr></thead><tbody>%s</tbody></table>', crows)
  # class totals table (Mha)
  trows<-paste(apply(ct,1,function(r){
    a<-km2_to_Mha(as.numeric(r[["actual"]])); m<-km2_to_Mha(as.numeric(r[["est_median"]]))
    lo<-km2_to_Mha(as.numeric(r[["est_q025"]])); hi<-km2_to_Mha(as.numeric(r[["est_q975"]]))
    pd<-as.numeric(r[["pdiff"]])
    sprintf('<tr><td class="cl">%s</td><td class="n">%.3f</td><td class="n">%.3f</td><td class="n ci">[%.3f, %.3f]</td><td class="n %s">%+.1f%%</td></tr>',
      esc(r[["class"]]), a, m, lo, hi, if(abs(pd)<2) "ok" else "warn2", pd)
  }),collapse="")
  tot<-sprintf('<table class="dt"><thead><tr><th>class</th><th class="n">actual (Mha)</th><th class="n">estimated (Mha)</th><th class="n">95%% CI</th><th class="n">Δ</th></tr></thead><tbody>%s</tbody></table>', trows)
  # per-country heat-table (estimated median Mha, actual on hover)
  cc<-as.data.table(R$ctry_tot); cc[,estM:=km2_to_Mha(est_median)]; cc[,actM:=km2_to_Mha(actual)]
  countries<-sort(unique(cc$group)); clsK<-cats[order(-colMeans(R$obs))]
  wide_est<-dcast(cc, group~class, value.var="estM"); setcolorder(wide_est, c("group",clsK))
  wide_act<-dcast(cc, group~class, value.var="actM"); setcolorder(wide_act, c("group",clsK))
  mx<-max(as.matrix(wide_est[,-1]),na.rm=TRUE)
  hh<-paste0("<th>region</th>", paste(sprintf('<th class="rot"><span>%s</span></th>',esc(clsK)),collapse=""))
  brows<-paste(sapply(seq_len(nrow(wide_est)),function(i){
    cells<-paste(sapply(clsK,function(cl){ v<-wide_est[[cl]][i]; av<-wide_act[[cl]][i]
      al<-if(is.na(v)||mx<=0) 0 else min(1,v/mx)
      sprintf('<td class="hm" style="background:rgba(15,122,103,%.3f)" title="%s / %s: est %.3f Mha (actual %.3f)">%s</td>',
        0.08+0.85*al, esc(wide_est$group[i]), esc(cl), if(is.na(v))0 else v, if(is.na(av))0 else av,
        if(is.na(v)) "" else sprintf("%.2f", v)) }),collapse="")
    sprintf('<tr><td class="cl sticky">%s</td>%s</tr>', esc(wide_est$group[i]), cells)
  }),collapse="")
  ctytab<-sprintf('<div class="scrollx"><table class="dt heat"><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>', hh, brows)
  # maps
  files<-list.files(file.path(MAPDIR,b), pattern="\\.png$", full.names=TRUE)
  dom<-files[grepl("00_dominant",files)]; heat<-files[grepl("heatplot",files)]
  clsfiles<-sort(files[grepl("^[0-9]",basename(files)) & !grepl("00_dominant",basename(files))])
  mapcards<-paste(sapply(clsfiles,function(f) sprintf('<figure class="mapc"><img loading="lazy" src="%s" alt=""></figure>', b64(f))),collapse="")
  header_block <- if (isTRUE(R$use_true_int)) {
    '<h3>Estimated vs actual class totals <span class="sub">genuine posterior</span></h3>
    <p class="note" style="border-left-color:var(--good)">Totals are posterior estimates from the recovered coefficients — <strong>exact slopes and intercepts</strong> (clean re-fit batches). The 95% intervals reflect the full posterior across draws.</p>'
  } else {
    '<h3>Estimated vs actual class totals <span class="cav">conditional reconstruction — see note</span></h3>
    <p class="note">Slopes are exact posterior draws; intercepts are re-solved conditional on them (the saved batches’ intercepts are corrupted). Because conditional intercepts absorb each region’s marginal composition, <strong>estimated totals track actual by construction</strong> and the 95% intervals here reflect slope-driven spatial uncertainty only (narrower than a full re-fit). Driver effects, spatial patterns and fit skill below are exact.</p>'
  }
  sprintf('
  <section class="branch" id="br-%s" %s>
    <div class="tiles">%s</div>
    %s
    <div class="two"><div>%s</div><div>%s</div></div>
    <h3>Totals per RE region <span class="sub">estimated median, Mha · hover for actual</span></h3>
    %s
    <h3>Driver effects <span class="sub">exact recovered slopes</span></h3>
    <figure class="wide"><img loading="lazy" src="%s" alt="driver effect heatplot"></figure>
    <h3>Spatial fit — dominant class</h3>
    <figure class="wide"><img loading="lazy" src="%s" alt="dominant class map"></figure>
    <h3>Observed vs predicted share — all %d classes</h3>
    <div class="mapgrid">%s</div>
  </section>', b, if(b==names(D)[1]) "" else "hidden", tiles, header_block,
    sprintf('<h3>Calibration <span class="sub">per-class share &amp; spatial correlation</span></h3>%s', cal),
    sprintf('<h3>Class totals <span class="sub">Mha</span></h3>%s', tot),
    ctytab, b64(heat), b64(dom), J, mapcards)
}
tabs<-paste(sprintf('<button class="tab%s" data-b="%s">%s<span class="tscheme">%s</span></button>',
  ifelse(names(D)==names(D)[1]," active",""), names(D), names(D), sapply(D,function(x)esc(x$cfg$scheme))), collapse="")
sections<-paste(sapply(names(D), branch_section), collapse="\n")

css <- '
:root{--bg:#FAFBFA;--surf:#fff;--ink:#14211E;--mut:#5E6E69;--line:#E5EAE7;--acc:#0F7A67;--acc2:#075447;--warn:#B26B00;--warnbg:#FBF3E4;--bad:#B2182B;--good:#0F7A67}
@media(prefers-color-scheme:dark){:root{--bg:#0E1613;--surf:#14201C;--ink:#E8F0EC;--mut:#93A39D;--line:#25332E;--acc:#54BCA6;--acc2:#8AD3C3;--warn:#E0A050;--warnbg:#241d10;--bad:#E27D8C}}
:root[data-theme=dark]{--bg:#0E1613;--surf:#14201C;--ink:#E8F0EC;--mut:#93A39D;--line:#25332E;--acc:#54BCA6;--acc2:#8AD3C3;--warn:#E0A050;--warnbg:#241d10;--bad:#E27D8C}
:root[data-theme=light]{--bg:#FAFBFA;--surf:#fff;--ink:#14211E;--mut:#5E6E69;--line:#E5EAE7;--acc:#0F7A67;--acc2:#075447;--warn:#B26B00;--warnbg:#FBF3E4;--bad:#B2182B;--good:#0F7A67}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;line-height:1.5;font-size:15px}
.wrap{max-width:1160px;margin:0 auto;padding:0 22px 80px}
header.top{position:sticky;top:0;z-index:20;background:color-mix(in srgb,var(--bg) 88%,transparent);backdrop-filter:blur(8px);border-bottom:1px solid var(--line)}
.hd{max-width:1160px;margin:0 auto;padding:14px 22px 0}
h1{font-size:22px;font-weight:640;letter-spacing:-.01em;margin:0 0 2px;text-wrap:balance}
.subttl{color:var(--mut);font-size:13.5px;margin:0 0 12px}
.tabs{display:flex;gap:4px;flex-wrap:wrap}
.tab{appearance:none;border:1px solid var(--line);border-bottom:none;background:transparent;color:var(--mut);padding:9px 15px 10px;border-radius:9px 9px 0 0;cursor:pointer;font:inherit;font-weight:600;display:flex;flex-direction:column;line-height:1.15}
.tab .tscheme{font-weight:400;font-size:10.5px;letter-spacing:.02em;opacity:.8}
.tab.active{color:var(--ink);background:var(--surf);border-color:var(--line);box-shadow:0 -2px 0 var(--acc) inset}
.intro{margin:26px 0 8px}
.caveat{border:1px solid color-mix(in srgb,var(--warn) 40%,var(--line));background:var(--warnbg);border-radius:12px;padding:14px 16px;margin:16px 0 8px;font-size:13.7px}
.caveat b{color:var(--warn)}
.branch{padding-top:20px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:6px 0 20px}
.tile{background:var(--surf);border:1px solid var(--line);border-radius:12px;padding:13px 14px}
.tval{font-size:23px;font-weight:650;letter-spacing:-.02em;font-variant-numeric:tabular-nums}
.tlabel{color:var(--mut);font-size:12px;margin-top:2px}
.tsub{color:var(--mut);font-size:11px;opacity:.8;margin-top:3px}
h3{font-size:15.5px;font-weight:620;margin:30px 0 10px;letter-spacing:-.005em}
h3 .sub,h3 .cav{font-weight:400;font-size:12px;color:var(--mut);margin-left:8px}
h3 .cav{color:var(--warn)}
.note{color:var(--mut);font-size:12.7px;background:var(--surf);border:1px solid var(--line);border-left:3px solid var(--warn);border-radius:8px;padding:10px 13px;margin:0 0 12px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:820px){.two{grid-template-columns:1fr}}
table.dt{width:100%;border-collapse:collapse;font-size:12.8px}
table.dt th,table.dt td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:left}
table.dt th{color:var(--mut);font-weight:600;font-size:11px;text-transform:uppercase;letter-spacing:.03em}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
td.cl{font-weight:500}
td.ci{color:var(--mut);font-size:11.5px}
td.warn2{color:var(--warn)} td.ok{color:var(--mut)}
.barcell{width:110px}.bar{display:inline-block;height:8px;border-radius:4px;background:var(--acc)}
.scrollx{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:#fff}
table.heat{font-size:11px;border-collapse:separate;border-spacing:0;background:#fff;color:#14211E}
table.heat th.rot{height:96px;white-space:nowrap;vertical-align:bottom;padding:0}
table.heat th.rot span{display:inline-block;transform:rotate(-90deg);transform-origin:left;translate:12px -6px;font-size:10px;color:#5E6E69}
table.heat td.hm{text-align:right;font-variant-numeric:tabular-nums;padding:3px 6px;color:#0a1512;border-bottom:1px solid rgba(255,255,255,.5)}
table.heat td.sticky{position:sticky;left:0;background:#fff;color:#14211E;font-weight:600}
figure{margin:0}
figure.wide img,figure.mapc img{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:10px;background:var(--surf)}
.mapgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:12px}
footer{color:var(--mut);font-size:12px;margin-top:40px;border-top:1px solid var(--line);padding-top:16px}
.themebtn{position:fixed;right:14px;bottom:14px;z-index:30;border:1px solid var(--line);background:var(--surf);color:var(--ink);border-radius:999px;padding:8px 12px;cursor:pointer;font:inherit;font-size:12px}
'
js <- "
document.querySelectorAll('.tab').forEach(t=>t.addEventListener('click',()=>{
 document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));t.classList.add('active');
 document.querySelectorAll('.branch').forEach(s=>s.hidden=true);
 document.getElementById('br-'+t.dataset.b).hidden=false;window.scrollTo({top:0,behavior:'smooth'});}));
const tb=document.getElementById('themebtn');function cur(){return document.documentElement.getAttribute('data-theme')|| (matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light')}
tb.addEventListener('click',()=>{const n=cur()==='dark'?'light':'dark';document.documentElement.setAttribute('data-theme',n);tb.textContent=n==='dark'?'☀ light':'☾ dark'});
"
n_cond <- sum(!sapply(D, function(x) isTRUE(x$use_true_int)))
caveat_banner <- if (n_cond > 0) sprintf('<div class="caveat"><b>Data-integrity note.</b> %s of the %d branch batch sets were written before a July-21 back-transform fix and carry a <b>corrupted intercept</b> (verified; the reconstructed log-likelihood misses the stored value entirely on the intercept alone). The regression <b>slopes recover exactly</b>, so driver effects, spatial patterns and fit skill shown here are trustworthy. For those branches intercepts were re-solved conditional on the exact slopes; <b>genuine class-total credible intervals require a re-fit</b> with the corrected sampler (verified clean, round-trip d≈0.2). Branches marked “genuine posterior” use clean re-fit batches.</div>',
  n_cond, length(D)) else '<div class="caveat" style="border-color:color-mix(in srgb,var(--good) 45%,var(--line))"><b>All branches from clean re-fit batches.</b> Slopes and intercepts are exact recovered posterior draws; class totals and credible intervals are genuine posterior estimates.</div>'
html <- sprintf('<title>Prior land-use model — fit report</title>
<style>%s</style>
<header class="top"><div class="hd"><h1>Prior land-use models — model-fitting report</h1>
<p class="subttl">10 km EU grid · zero-sum Bayesian multinomial (RE by region) · AGMIP / GLOBIOM / BIOCLIMA branches</p>
<div class="tabs">%s</div></div></header>
<div class="wrap">
%s
%s
<footer>Reconstructed from on-disk posterior batches (slopes exact, %s draws/branch) · conditional-intercept prediction · generated %s. Totals in Mha (1 Mha = 10,000 km²).</footer>
</div>
<button class="themebtn" id="themebtn">☾ dark</button>
<script>%s</script>', css, tabs, caveat_banner, sections, D[[1]]$ndraw_used, format(Sys.Date()), js)

writeLines(html, OUTHTML)
cat("wrote", OUTHTML, "\n"); cat("size:", round(file.info(OUTHTML)$size/1e6,1), "MB\n")
