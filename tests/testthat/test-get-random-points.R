test_that ('generate_random_points', {
               expect_error (generate_random_points (),
                             "city is missing with no default.")
               expect_error (generate_random_points (city = "Barcelona"),
                             "n is missing with no default.")
               expect_error (generate_random_points (n = 100),
                             "city is missing with no default.")
})
