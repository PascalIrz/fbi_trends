FBI trends
================
2024-03-28

# Aim

This repo contains the code to analyse the temporal trends in ‘river
health’ using a Fish-based Index (FBI; see Oberdorff *et al.* 2001,
2002). This is intended to submit an article.

# Repo structure

Main sub-directories are:

1.  `scripts`: contains a series of R and R Markdown scripts to be
    executed in the order they are numbered. The complete workflow can
    be done using the file *99_make.R*. Before running the scripts, the
    database is required and can be found
    [here](https://doi.org/10.5281/zenodo.7099129) (See Irz *et al.*
    2022). Note that the `refnet_start` parameter (first year of the
    REFNET data set) that is passed to some of the R Markdown files.
2.  `R`: contains some custom functions for the analysis of the data and
    the graphical outputs.

The file `60_output.html` contains the elements for the article.

# Main references

- Oberdorff T, Pont D, Hugueny B, Chessel D. 2001. A probabilistic model
  characterizing fish assemblages of French rivers: A framework for
  environmental assessment: Predicting riverine fish assemblages.
  *Freshw Biol* 46: 399–415.

- Oberdorff T, Pont D, Hugueny B, Porcher J-P. 2002. Development and
  validation of a fish-based index for the assessment of “river health”
  in France. *Freshw Biol* 47: 1720–1734.

- Irz P, Vigneron T, Poulet N, Cosson E, Point T, Baglinière E, Porcher
  JP. 2022. A long-term monitoring database on fish and crayfish species
  in French rivers. *Knowl. Manag. Aquat. Ecosyst.* 423, 25.
