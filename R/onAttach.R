SPLICE_EVENTS = list(
  "SkippingExon" = "SE",
  "Alternative5Prime" = "A5",
  "Alternative3Prime" = "A3",
  "MutuallyExclusiveExon" = "MX",
  "RetainedIntron" = "RI",
  "AlternativeFirstExon" = "AF",
  "AlternativeLastExon" = "AL",
  "Isoform" = 'isoform'
)

SPLICE_EVENTS.DECODE = unlist(SPLICE_EVENTS)
SPLICE_EVENTS.REVERSE = names(SPLICE_EVENTS)
names(SPLICE_EVENTS.REVERSE) = unlist(SPLICE_EVENTS)

SUPPA_PATH = "/slipstream/home/joeboyd/anaconda2/envs/suppa2_env/bin/suppa.py"



.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching ssvSplicing version ",
                        packageDescription("ssvSplicing")$Version, ".")
  SSV_SPLICE_EVENTS <<- SPLICE_EVENTS
  SSV_SPLICE_EVENTS.DECODE <<- SPLICE_EVENTS.DECODE
  SSV_SPLICE_EVENTS.REVERSE <<- SPLICE_EVENTS.REVERSE
}
