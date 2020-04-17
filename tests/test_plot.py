#!/usr/bin/env python3
"""PyTest tests for the plot.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from plot import *



def test_get_text_label_matrix():
    """Test the get_text_label_matrix function in the plot.py file.
    """
    ##########################
    # Arrange.
    odf = 'odf'
    sequence_type = 'sequence_type'

    ##########################
    # Act.
    #x = get_text_label_matrix(odf,
    #		sequence_type)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_count_matrix():
    """Test the get_hit_count_matrix function in the plot.py file.
    """
    ##########################
    # Arrange.
    odf = 'odf'
    sequence_type = 'sequence_type'

    ##########################
    # Act.
    #x = get_hit_count_matrix(odf,
    #		sequence_type)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_heat_matrix_from_count_matrix():
    """Test the get_heat_matrix_from_count_matrix function in the plot.py file.
    """
    ##########################
    # Arrange.
    data_count = 'data_count'

    ##########################
    # Act.
    #x = get_heat_matrix_from_count_matrix(data_count)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_generate_heatmap_figure():
    """Test the generate_heatmap_figure function in the plot.py file.
    """
    ##########################
    # Arrange.
    column_labels = 'column_labels'
    row_labels = 'row_labels'
    data_labels = 'data_labels'
    data_count = 'data_count'

    ##########################
    # Act.
    #x = generate_heatmap_figure(column_labels,
    #		row_labels,
    #		data_labels,
    #		data_count)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_generate_stats_file():
    """Test the generate_stats_file function in the plot.py file.
    """
    ##########################
    # Arrange.
    column_labels = 'column_labels'
    row_labels = 'row_labels'
    data_count = 'data_count'

    ##########################
    # Act.
    #x = generate_stats_file(column_labels,
    #		row_labels,
    #		data_count)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_sector_label():
    """Test the get_sector_label function in the plot.py file.
    """
    ##########################
    # Arrange.
    pct = 'pct'
    iterable = 'iterable'

    ##########################
    # Act.
    #x = get_sector_label(pct,
    #		iterable)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_complex_info_list():
    """Test the get_complex_info_list function in the plot.py file.
    """
    ##########################
    # Arrange.
    complex_info_file = 'complex_info_file'
    column_labels_simple = 'column_labels_simple'

    ##########################
    # Act.
    #x = get_complex_info_list(complex_info_file,
    #		column_labels_simple)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_modify_legend_lines():
    """Test the modify_legend_lines function in the plot.py file.
    """
    ##########################
    # Arrange.
    wedges = 'wedges'

    ##########################
    # Act.
    #x = modify_legend_lines(wedges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_modify_lines():
    """Test the modify_lines function in the plot.py file.
    """
    ##########################
    # Arrange.
    wedges = 'wedges'

    ##########################
    # Act.
    #x = modify_lines(wedges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_make_coulson_plot():
    """Test the make_coulson_plot function in the plot.py file.
    """
    ##########################
    # Arrange.
    column_labels_simple = 'column_labels_simple'

    ##########################
    # Act.
    #x = make_coulson_plot(column_labels_simple)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_plot_amoebae_res():
    """Test the plot_amoebae_res function in the plot.py file.
    """
    ##########################
    # Arrange.
    csv_file = 'csv_file'
    complex_info = 'complex_info'
    outpdfpath = 'outpdfpath'
    csv_file2 = 'csv_file2=None'

    ##########################
    # Act.
    #x = plot_amoebae_res(csv_file,
    #		complex_info,
    #		outpdfpath,
    #		csv_file2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


