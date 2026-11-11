# Shared helpers used by the final fixed-study analysis.

outcome_sampled_reinclusion_rows <- function() {
  rbind(
    make_arm_rows(
      "Iachetta_2019",
      c("13hetho", "2846hetho", "WT_Biased"),
      c(2, 5, 361), c(2, 5, 141)
    ),
    make_arm_rows(
      "Kleibl_2009",
      c("HapB3_1129or1236", "2Ahetho", "WT_Biased"),
      c(4, 5, 119), c(2, 5, 71)
    )
  )
}
