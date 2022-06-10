--- 
title: "Learning to Love dplR"
subtitle: "Or Using R for Tree-Ring Analysis"
author: "Andy Bunn"
date: "10-June-2022"
site: bookdown::bookdown_site
output: 
  bookdown::gitbook:
    config:
      toc:
        collapse: section
        scroll_highlight: true
        before: null
        after: null
    split_bib: no
documentclass: book
bibliography: [dplR.bib, packages.bib]
csl: global-change-biology.csl
biblio-style: apalike
link-citations: yes
github-repo: OpenDendro/dplR-workshop
description: "Helful materials for dplR"
---


# Preamble

`R` is both a programming language and a software environment for statistical computing. It is free and open-source, licensed under the GNU General Public License. The Dendrochronology Program Library in R (`dplR`) is an add-on package for R that performs many of the standard tasks in tree-ring analysis including cross-dating, detrending, chronology building, spectral and wavelet analyses, and so on.  

In this workshop, we will use example ring-width files from the ITRDB to demonstrate the functionality of dplR. The R environment is powerful and flexible and its use allows great transparency in presenting data results. An advantage of dplR's open-source licensing is that one has the option to modify or add to the library's functionality for performing specific experiments or producing custom figures.

In the following pages we will cover the basics of `dplR`, crossdating, and some limited time-series analysis. Users should be familiar with the basics of dendrochronology and concepts like detrending, autocorrelation, spectral analysis and so on.

If this is all new to you -- you should proceed immediately to a good primer on dendrochronology like Fritts [-@Fritts2001] or the Cook Book [-@CookBook]. These pages are not intended to teach you about how to do tree-ring analysis. They are intended to teach you how to use R for dendro.

**Please note! This is a very drafty document and there are typos and all kinds of silliness in it.**

This document was written in **Markdown** using the  `bookdown` package.

