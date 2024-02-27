def scale_bar_image(image, 
                    micron_length, 
                    scale_length,
                    scale_color="white",
                    thickness=4,
                    bottom_margin=0.93,
                    left_margin=0.07):                   
    
    """
    Generates a image with a scale bar at the bottom left side
    
    Parameters:
    - image (numpy.ndarray): The input image (2D numpy array) which will be used to generate the scale bar image.
    - micron_length (float): The physical length (in microns) represented by the scale bar.
    - scale_length (float): The desired length of the scale bar in the generated image.
    - scale_color (str, optional): Color of the scale bar. Default is "white".
    - thickness (int, optional): Thickness of the scale bar. Default is 4 pixels.
    - bottom_margin (float, optional): The position of the scale bar from the bottom as a fraction of the image height. Default is 0.93.
    - left_margin (float, optional): The position of the scale bar from the left as a fraction of the image width. Default is 0.07.
    
    Returns:
    - scale_bar (numpy.ndarray): The image with the overlaid scale bar.
    - colormap (matplotlib.colors.ListedColormap): Colormap used for the scale bar.

    Raises:
    - ValueError: If the input image does not have 2 dimensions.
    
    Example:
    
    #import clesperanto
    import pyclesperanto_prototype as cle
        
    #load image
    image = cle.imread("file_path)
    
    #create scale bar and colormap
    scale_bar, colormap = scale_bar_image(image, micron_length=0.323, scale_length=25)
    
    #show image
    cle.imshow(image, continue_drawing=True)
    
    #overlay scale bar
    cle.imshow(scale_bar, colormap=colormap)
    """
    
    if len(image.shape) != 2:
        raise ValueError("Out of index. 2-dimensional images are required.")
    
    from matplotlib.colors import ListedColormap
    import numpy as np
    
    colormap = ListedColormap(["none", scale_color])
    
    scale_bar = np.copy(image)
    scale_bar[:] = 0
    
    length_pixels = int(scale_length / micron_length)
    
    margin_bottom = int(scale_bar.shape[0] * bottom_margin)
    margin_left = int(scale_bar.shape[1] * left_margin)
    
    scale_bar[margin_bottom:margin_bottom + thickness, 
              margin_left:margin_left + length_pixels] = 1
    
    return scale_bar, colormap

def correlation_matrix(dataframe, method="pearson", font_scale=0.75):
    """
    Generate a heatmap representing the correlation matrix of a given dataframe.

    Parameters:
    - dataframe (pd.DataFrame): The input dataframe containing numerical data.
    - method (str, optional): The method used to compute correlation. Default is "pearson".
    - font_scale: The font_size of the correlation matrix. Default is 0.75.

    Returns:
    - seaborn.heatmap: A seaborn heatmap representing the correlation matrix.

    Dependencies:
    - seaborn
    - matplotlib.colors.LinearSegmentedColormap
    - numpy (as np)

    Example:
    >>> import pandas as pd
    >>> import numpy as np
    >>> import seaborn as sns
    >>> df = pd.DataFrame(np.random.randn(10, 5), columns=['A', 'B', 'C', 'D', 'E'])
    >>> correlation_matrix(df)
    """
    import numpy as np
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    
    # compute the correlation matrix and respective mask
    correlation_matrix = dataframe.corr(method=method)
    mask = np.triu(correlation_matrix)
    
    # compute a color gradient for the correlation matrix
    colors = [(0, 0, 0.5), (1, 1, 1), (0.5, 0, 0)]
    gradient = LinearSegmentedColormap.from_list('custom_gradient', colors, N=256)
    
    # create the plot
    sns.set(font_scale=font_scale)
    heatmap = sns.heatmap(data=correlation_matrix,
                          mask=mask,
                          center=0.0,
                          vmin=-1.0,
                          vmax=1.0,
                          cmap=gradient,
                          square=True)
    
    return heatmap

def correlation_vector_age(dataframe, method="pearson", font_scale=0.75):
    
    """
    Compute the correlation vector of age in the given dataframe and visualize it as a heatmap.

    Parameters:
    - dataframe (pd.DataFrame): The pandas DataFrame containing the data.
    - method (str, optional): The method to compute correlation. Default is "pearson".
    - font_scale (float, optional): The font size of the correlation vector. Default is 0.75.

    Returns:
    - seaborn.matrix.ClusterGrid: Heatmap visualization of the correlation vector of age.
    """
    
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    
    #compute correlation vector of the age
    correlation_matrix = dataframe.corr(method=method)
    correlation_vector_age = correlation_matrix[["age"]]
    
    #sort the values of the correlation vector
    correlation_vector_age = correlation_vector_age.sort_values(by="age", ascending=False)
    correlation_vector_age = correlation_vector_age.drop(labels="age", axis=0)
    correlation_vector_age = correlation_vector_age.round(2)
    
    #define a color gradient
    colors = [(0, 0, 0.5), (1, 1, 1), (0.5, 0, 0)]
    gradient = LinearSegmentedColormap.from_list('custom_gradient', colors, N=256)
    
    #compute the heatmap
    sns.set(font_scale = font_scale)
    heatmap = sns.heatmap(data=correlation_vector_age,
                          cmap=gradient,
                          vmin=-1.0,
                          vmax=1.0,
                          center=0.0,
                          annot=True,
                          square=True)
    
    return heatmap

def significant_tukey_results(dataframe, x, y):
    """
    Computes Tukey's Honestly Significant Difference (HSD) test on the provided dataframe.
    
    Args:
        dataframe (pd.DataFrame): The dataframe containing the data for analysis.
        x (str): The name of the column representing the independent variable.
        y (str): The name of the column representing the dependent variable.
        
    Returns:
        pd.DataFrame: A dataframe containing significant Tukey results sorted by the second group.
        
    Note:
        Before using this function, ensure that the column 'y' has been averaged by the column 'x'.
        This function is typically used to annotate significant differences in groups on a line plot.
    """
    import pandas as pd
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    
    # Compute Tukey Results
    tukey_results = pairwise_tukeyhsd(dataframe[y], dataframe[x])
    
    # Turn tukey results into a pandas dataframe
    tukey_results = pd.DataFrame(data=tukey_results._results_table.data[1:], 
                                  columns=tukey_results._results_table.data[0])
    
    # Retrieve the significant rows only and sort them by the second group of significance
    significant_rows = tukey_results[tukey_results["reject"] == True]
    significant_rows = significant_rows.sort_values(by="group2")
    
    return significant_rows

def x_coordinate(xmin, xmax):
    """
    Determine the x-coordinate of a significance star in a plot.

    Parameters:
    - xmin (float): The minimum x-value of the plot.
    - xmax (float): The maximum x-value of the plot.

    Returns:
    float: The x-coordinate of the significance star, which is located halfway between xmin and xmax.
    """
    # Determine the x position
    x_coord = xmin + ((xmax - xmin) / 2)
    
    return x_coord

def asterics(row, column_name):
    """
    Convert p-values into significance stars.

    Parameters:
    - row (pandas.Series): A row from a DataFrame containing at least a column named "p-adj" representing p-values.

    Returns:
    - str: A string representing the significance stars based on the p-value.
      - "***" if p-value < 0.001
      - "**" if 0.001 <= p-value < 0.01
      - "*" if 0.01 <= p-value < 0.05
    """
    # Determine p value
    p_value = row[column_name]
    
    # Determine the significance star
    if p_value < 0.001:
        significant_result = "***"
    elif p_value < 0.01:
        significant_result = "**"
    elif p_value < 0.05:
        significant_result = "*"
    
    return significant_result

def x_coordinate_barplot(row, ages):
    """
    Determine the x-coordinate for drawing an asterisk on a barplot based on the provided row and ages array.

    Parameters:
    - row (dict): A dictionary representing a row from a t_test_table containing at least the "age" key.
    - ages (list): A list of ages used for comparison.

    Returns:
    - int: The x-coordinate where the asterisk should be drawn on the barplot.
    """
    
    # Iterate through the ages to find the matching age in the row
    for x_coordinate, age in enumerate(ages):
        
        # If the age in the row matches the current age in the iteration, return the x-coordinate
        if age == row["age"]:
            return x_coordinate

def significant_t_tests(dataframe, ages, y, x="age", wt="wt", test="cpfl", hue="genotype"):
    """
    Perform significant t-tests between two groups at different ages.

    Parameters:
    -----------
    dataframe : pandas DataFrame
        The DataFrame containing the data.
    ages : list-like
        A list or array-like object containing the different ages at which tests are conducted.
    y : str
        The column name in the DataFrame representing the measurement of interest.
    x : str, optional
        The column name in the DataFrame representing the age variable (default is 'age').
    wt : str, optional
        The label representing the wild-type group in the 'hue' column (default is 'wt').
    test : str, optional
        The label representing the test sample group in the 'hue' column (default is 'cpfl').
    hue : str, optional
        The column name in the DataFrame representing the categorical variable distinguishing groups (default is 'genotype').

    Returns:
    --------
    pandas DataFrame
        A DataFrame containing the significant t-test results, with columns 'age' and 'p_value'.

    Notes:
    ------
    - This function performs independent t-tests between the 'wt' and 'test' groups for each specified age.
    - It filters the results based on the significance level of the t-tests (p < 0.05).
    """
    
    # Import necessary packages
    import pandas as pd
    from scipy.stats import ttest_ind
    
    # Define empty array to store test results
    test_results = []
    
    # Test different ages with independent t test
    for postnatal_age in ages:
        
        # Filter a specific age
        age_filtered = dataframe[dataframe[x] == postnatal_age]
    
        # Filter the measurement of interest in wt and test sample
        measurement_test_y = age_filtered[age_filtered[hue] == test][y]
        measurement_wt_y = age_filtered[age_filtered[hue] == wt][y]
        
        # Perform the t test of the 2 groups
        t_statistic, p_value = ttest_ind(measurement_wt_y,
                                         measurement_test_y,
                                         equal_var=False,
                                         alternative="greater")
        
        # Append the test result into the empty dataframe, two columns age and p_value
        test_results.append({"age": postnatal_age, "p_value": p_value})
    
    # Turn into pandas dataframe and filter significant result
    test_results = pd.DataFrame(test_results)
    test_results = test_results[test_results["p_value"] < 0.05]
    
    return test_results

def normalise_image(image):
    '''
    Normalizes an intensity image for training data used by the OS pixel classifier (quapos-lm.cl).

    Parameters:
    image (numpy.ndarray): Input intensity image to be normalized.

    Returns:
    numpy.ndarray: Normalized intensity image.

    Normalizes the input intensity image by performing the following steps:
    1. Applies top-hat background subtraction using a box-shaped structuring element with a radius of 3 pixels.
    2. Calculates the maximum intensity value in the background-subtracted image.
    3. Normalizes the intensity values by dividing by the maximum intensity value and scaling to fit within the dynamic range [0, 2^12].
    
    Example:
    >>> normalized_image = normalise_image(intensity_image)
    '''
    
    import numpy as np
    import pyclesperanto_prototype as cle
    
    # Apply top-hat background subtraction to the input image
    image_background_subtracted = cle.top_hat_box(source=image,
                                                  radius_x=3,
                                                  radius_y=3,
                                                  radius_z=3)
    
    # Find the maximum intensity value in the background-subtracted image
    max_intensity = np.max(image_background_subtracted)
    
    # Normalize the intensity values to fit within [0, 2^12] dynamic range
    image_normalised_intensity = (image_background_subtracted / max_intensity) * np.power(2, 12)
    
    return image_normalised_intensity

def heatmap_image(segmentation, dataframe, feature):
    """
    Generate a heatmap image based on segmentation and a specified feature from a dataframe for annotation.

    Parameters:
    - segmentation: Segmentation array representing regions of interest.
    - dataframe: DataFrame containing data for annotation.
    - feature (str): Name of the feature/column in the dataframe to be used for annotation.

    Returns:
    - pyclesperanto_prototype._tier0._pycl.OCLArray: Heatmap image generated by replacing intensities in the segmentation based on the specified feature.
    
    """
    # Import packages
    import pyclesperanto_prototype as cle
    
    # Retrieve a list from the feature for annotation
    feature_list = [0] + dataframe[feature].tolist()
    
    # Turn prediction into a heatmap image
    heatmap_image = cle.replace_intensities(segmentation, feature_list)
    
    # Return the heatmap image
    return heatmap_image

def predict_image(image, classifier):
    """
    Predict the labels of a microscopy image using a given machine learning classifier.

    Args:
        image (numpy.ndarray): The microscopy image to be analyzed.
        classifier (apoc.model): The machine learning classifier used for prediction. 
            Defaults to quapos_lm.

    Returns:
        list: The predicted labels for the input image.

    Raises:
        ImportError: If the required modules 'apoc' or 'fluorescent_microscopy_analysis' 
            are not installed or cannot be imported.
    """
    import apoc
    from quapos_lm import normalise_image

    # Normalize the input image using the provided function
    normalised_image = normalise_image(image)

    # Predict labels for the normalized image using the provided classifier
    prediction = classifier.predict(normalised_image)

    return prediction

def rescale_image(image, voxel_x, voxel_y, voxel_z):
    """
    Rescale an anisotropic image to isotropic voxel dimensions.

    Parameters:
    -----------
    image : cle.Image
        The input anisotropic image to be rescaled.
    voxel_x : float
        Voxel size along the x-axis in the original image (in units).
    voxel_y : float
        Voxel size along the y-axis in the original image (in units).
    voxel_z : float
        Voxel size along the z-axis in the original image (in units).

    Returns:
    --------
    cle.Image
        The rescaled image where all voxel sizes are equal.

    Notes:
    ------
    This function rescales an anisotropic image to make all voxel dimensions equal.
    It uses pyclesperanto_prototype.scale() function for rescaling.

    Example:
    --------
    >>> rescaled_image = rescale_image(anisotropic_image, 0.5, 0.8, 1.2)
    """
    
    import pyclesperanto_prototype as cle
    
    # Calculate zoom factor based on the smallest voxel dimension
    zoom_factor = 1 / voxel_x
    
    # Scale factors for each axis
    scale_x = voxel_x * zoom_factor
    scale_y = voxel_y * zoom_factor
    scale_z = voxel_z * zoom_factor
    
    # Rescale the image using pyclesperanto_prototype.scale()
    image_rescaled = cle.scale(source=image,
                               factor_x=scale_x,
                               factor_y=scale_y,
                               factor_z=scale_z,
                               auto_size=True)
    
    return image_rescaled

def rescale_segmentation(segmentation, voxel_x, voxel_y, voxel_z):
    """
    Rescale an anisotropic segmentation to isotropic voxel dimensions.

    Parameters:
    -----------
    segmentation: cle.Image
        The input anisotropic segmentation to be rescaled.
    voxel_x : float
        Voxel size along the x-axis in the original image (in units).
    voxel_y : float
        Voxel size along the y-axis in the original image (in units).
    voxel_z : float
        Voxel size along the z-axis in the original image (in units).

    Returns:
    --------
    cle.Image
        The rescaled segmentation where all voxel sizes are equal.

    Notes:
    ------
    This function rescales an anisotropic image to make all voxel dimensions equal.
    It uses pyclesperanto_prototype.scale() function for rescaling.

    Example:
    --------
    >>> rescaled_segmentation = rescale_segmentation(anisotropic_segmentation, 0.5, 0.8, 1.2)
    """
    
    import pyclesperanto_prototype as cle
    
    # Calculate zoom factor based on the smallest voxel dimension
    zoom_factor = 1 / voxel_x
    
    # Scale factors for each axis
    scale_x = voxel_x * zoom_factor
    scale_y = voxel_y * zoom_factor
    scale_z = voxel_z * zoom_factor
    
    # Rescale the segmentation using pyclesperanto_prototype.scale()
    segmentation_rescaled = cle.scale(source=segmentation,
                                      factor_x=scale_x,
                                      factor_y=scale_y,
                                      factor_z=scale_z,
                                      auto_size=True)
    
    # Cle.scale function changes the type of the image, change back into labelled result
    binary = segmentation_rescaled > 0
    
    # Relabel image
    segmentation_rescaled = cle.connected_components_labeling_diamond(binary)
    
    return segmentation_rescaled