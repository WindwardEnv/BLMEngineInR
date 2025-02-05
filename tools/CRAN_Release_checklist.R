devtools::load_all()

# update data files
sapply(list.files("data-raw", full.names = TRUE), FUN = source)

# check lints and style
styler::style_pkg(
  style = styler::tidyverse_style(scope = c("indention", "line_breaks")),
  dry = "on"
)
lintr::lint_package()

# to-do list before release:
# release()

# Update documentation
devtools::document()

# run tests and examples
devtools::run_examples()
devtools::test()

# check test coverage
devtools::unload()
tmp = devtools::test_coverage()
devtools::load_all()

# update manual
devtools::build_manual()

# Have you checked for spelling errors (with `spell_check()`)?
devtools::spell_check()
# update inst/WORDLIST with anything that might be flagged as misspelled or with
# spelling::update_wordlist()

# Have you run `R CMD check` locally?
devtools::check()
# make sure there are no errors, warnings, or notes...nothing...nada...get it?

# Were devtool's checks successful?
devtools::release_checks()
# These should all say "OK"

# Have you checked on R-hub (with `check_rhub()`)?
# check_rhub()
# This function is deprecated and defunct since rhub v2.
# Please see `?rhubv2` on transitioning to the new rhub functions.
rhub::rhub_check(platforms = c("linux", "macos", "macos-arm64", "windows"))

# Have you checked on win-builder (with `check_win_devel()`)?
check_win_devel()

# Have you updated `NEWS.md` file?
# Have you updated `DESCRIPTION`?
# Have you updated `cran-comments.md?`
# Were Git checks successful?
# Is your email address kellyc@windwardenv.com?
