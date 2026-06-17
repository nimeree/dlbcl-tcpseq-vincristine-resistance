# Run this script to capture the R session information for reproducibility.
# Output is saved to session_info.txt in the repository root.
sink("session_info.txt")
date()
sessionInfo()
sink()
message("Session info written to session_info.txt")
