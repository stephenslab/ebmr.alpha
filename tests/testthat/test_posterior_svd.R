test_that("posterior computations via woodbury and svd of Xtilde match direct ones",
{
            # testing
            n = 10
            p = 20
            X = matrix(rnorm(n * p), nrow = n)
            w = rnorm(p) ^ 2
            Xtilde = t(w ^ 0.5 * t(X))
            Xtilde.svd = svd(Xtilde)
            sb2 = 0.5

            expect_equal(Sigma1_direct(w, X, sb2),
                      Sigma1_woodbury_svd(w, Xtilde.svd, sb2))
            expect_equal(diag(Sigma1_direct(w, X, sb2)),
                      Sigma1_diag_woodbury_svd(w, Xtilde.svd, sb2))

            sb2 = 2

            expect_equal(Sigma1_direct(w, X, sb2),
                      Sigma1_woodbury_svd(w, Xtilde.svd, sb2))
            expect_equal(diag(Sigma1_direct(w, X, sb2)),
                      Sigma1_diag_woodbury_svd(w, Xtilde.svd, sb2))

            y = rnorm(n)
            expect_equal(mu1_direct(y, w, X, sb2), mu1_woodbury_svd(y, w, Xtilde.svd, sb2))
})
