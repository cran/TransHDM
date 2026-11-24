toy_frame <- function(n = 6) {
  data.frame(
    Y = seq_len(n) + 0.1,
    D = seq_len(n) / n,
    M1 = seq_len(n) * 0.3,
    M2 = rev(seq_len(n)) * 0.2,
    X1 = (seq_len(n) %% 3) + 0.4
  )
}

fake_transhdm <- function(n_med = 2L) {
  mediators <- paste0("M", seq_len(max(n_med, 1L)))
  contributions <- data.frame(
    mediator = mediators,
    alpha = c(0.5, -0.2)[seq_len(length(mediators))],
    alpha_pv = c(0.01, 0.2)[seq_len(length(mediators))],
    beta = c(0.4, 0.3)[seq_len(length(mediators))],
    beta_pv = c(0.02, 0.04)[seq_len(length(mediators))],
    alpha_beta = c(0.2, -0.06)[seq_len(length(mediators))],
    ab_pv = c(0.02, 0.2)[seq_len(length(mediators))],
    fdr_est = c(0.01, 0.04)[seq_len(length(mediators))],
    pa = c(0.1, -0.03)[seq_len(length(mediators))],
    stringsAsFactors = FALSE
  )
  if (n_med == 0L) {
    contributions <- contributions[0, , drop = FALSE]
  } else {
    contributions <- contributions[seq_len(n_med), , drop = FALSE]
  }

  effects <- data.frame(
    effect = c("indirect", "direct", "total", "pe"),
    estimate = c(0.14, 0.5, 0.64, 0.14 / 0.64),
    stringsAsFactors = FALSE
  )
  all_mediators <- data.frame(
    mediator = c("M1", "M2", "M9"),
    alpha = c(0.5, -0.2, 0.1),
    alpha_pv = c(0.01, 0.2, 0.4),
    beta = c(0.4, 0.3, 0),
    beta_pv = c(0.02, 0.04, 0.8),
    significant = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      contributions = contributions,
      effects = effects,
      IDE_est = c(M1 = 0.2, M2 = -0.06, M9 = 0),
      DE_est = c(D = 0.5),
      all_mediators = all_mediators
    ),
    class = "TransHDM"
  )
}

fake_inference <- function(mediators = c("M1", "M2")) {
  structure(
    list(
      DE_est = 0.4,
      alpha_SIS_est = setNames(c(0.5, 0.1)[seq_along(mediators)], mediators),
      beta_SIS_est = setNames(c(0.3, 0.2)[seq_along(mediators)], mediators),
      P_alpha_SIS = setNames(c(0.01, 0.4)[seq_along(mediators)], mediators),
      P_beta_SIS = setNames(c(0.02, 0.5)[seq_along(mediators)], mediators),
      M = mediators,
      M_SIS = mediators,
      D = "D"
    ),
    class = "mediation_inference"
  )
}

fake_sis <- function(with_source = FALSE) {
  target <- toy_frame()[, c("Y", "D", "M1", "X1")]
  structure(
    list(
      target_SIS = target,
      source_SIS = if (with_source) target else NULL,
      M_ID_name_SIS = "M1",
      Y = "Y",
      D = "D",
      M = c("M1", "M2"),
      X = "X1"
    ),
    class = "SIS"
  )
}

fake_source <- function(ids = 1L) {
  structure(
    list(
      transfer.source.id = ids,
      source.loss = c(1.2, 2.4),
      target.valid.loss = 1.5,
      T_index = c(-0.3, 0.9),
      threshold = 0.03,
      loss.cv = matrix(c(1.1, 1.3, 2.3, 2.5, 1.4, 1.6), nrow = 2)
    ),
    class = "source_detection"
  )
}
