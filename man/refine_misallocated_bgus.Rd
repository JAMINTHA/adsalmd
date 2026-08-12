\name{refine_misallocated_bgus}
\alias{refine_misallocated_bgus}
\title{Local search refinement: reassign misallocated BGUs to a better-fitting}
\usage{
refine_misallocated_bgus(best_sol, W, adj, row_W, col_W, max_iter = 100)
}
\arguments{
\item{best_sol}{Integer vector, starting partition to improve.}
\item{W}{Numeric matrix (N x N).}
\item{adj}{Integer/logical matrix (N x N).}
\item{row_W}{Numeric vector, rowSums(W).}
\item{col_W}{Numeric vector, colSums(W).}
\item{max_iter}{Integer, maximum number of full passes over all BGUs (default 100).}
}
\value{
Integer vector, the improved (relabelled) partition.
}
\description{
adjacent cluster. Repeatedly scans every BGU. If a BGU is misallocated (per is_bgu_misallocated), it is reassigned to whichever adjacent cluster gives the highest compute_ci score, as long as the resulting partition keeps every cluster contiguous. Passes repeat until a full pass makes no further improvement, or max_iter passes have been made.
}
