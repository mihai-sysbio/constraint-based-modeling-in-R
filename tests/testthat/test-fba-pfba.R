test_that("FBA optimizes toy model using GLPK backend", {
  model <- make_minimal_flux_model()

  result <- fba(model)

  expect_s4_class(result, "FluxPrediction")
  expect_identical(result@algorithm, "FBA")
  expect_identical(result@ok_term, "optimization process was successful")
  expect_identical(result@stat_term, "solution is optimal")
  expect_equal(result@obj, 10, tolerance = 1e-6)
  expect_equal(result@fluxes, c(10, 10), tolerance = 1e-6)
})

test_that("pFBA preserves objective and minimizes weighted flux", {
  model <- make_minimal_flux_model()

  baseline <- fba(model)
  result <- pfba(model)

  expect_s4_class(result, "FluxPrediction")
  expect_equal(result@obj, baseline@obj, tolerance = 1e-6)
  expect_equal(result@fluxes, baseline@fluxes, tolerance = 1e-6)
  expect_equal(result@obj_sec, sum(abs(baseline@fluxes)), tolerance = 1e-6)

  weighted <- pfba(
    model,
    costcoeffw = c(1, 3),
    costcoefbw = c(2, 4)
  )
  expect_equal(
    weighted@obj_sec,
    sum(baseline@fluxes * c(1, 3)),
    tolerance = 1e-6
  )

  expect_error(
    pfba(model, costcoeffw = c(1)),
    "must both be of length equal to the number of reactions"
  )
})
