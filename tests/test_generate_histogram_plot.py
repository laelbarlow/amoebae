#!/usr/bin/env python3
"""PyTest tests for the generate_histogram_plot.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from generate_histogram_plot import *



def test_generate_histogram():
    """Test the generate_histogram function in the generate_histogram_plot.py file.
    """
    ##########################
    # Arrange.
    title = 'title'

    ##########################
    # Act.
    #x = generate_histogram(title)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_generate_double_histogram():
    """Test the generate_double_histogram function in the generate_histogram_plot.py file.
    """
    ##########################
    # Arrange.
    title = 'title'

    ##########################
    # Act.
    #x = generate_double_histogram(title)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_autolabel_bars():
    """Test the autolabel_bars function in the generate_histogram_plot.py file.
    """
    ##########################
    # Arrange.
    rects = 'rects'
    ax = 'ax'

    ##########################
    # Act.
    #x = autolabel_bars(rects,
    #		ax)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_generate_bar_chart():
    """Test the generate_bar_chart function in the generate_histogram_plot.py file.
    """
    ##########################
    # Arrange.
    title = 'title'

    ##########################
    # Act.
    #x = generate_bar_chart(title)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


