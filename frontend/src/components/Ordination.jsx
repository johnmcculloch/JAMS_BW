import React, {useEffect, useState } from 'react';
import Button from '@mui/material/Button';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';

const Ordination = ({ handleNavigateTo }) => {
    // Parameter options
    const [parameters, setParameters] = useState({
        glomby: '',
        algorithm: '',
        PCA_Components: '',
        distmethod: '',
        subsetby: '',
        compareby: '',
        colourby: '',
        colorby: '',
        shapeby: '',
        samplesToKeep: '',
        featuresToKeep: '',
        samplesToHighlight: '',
        ignoreunclassified: true,
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        discard_SDoverMean_below: '',
        asPPM: true,
        normalization: 'relabund',
        PPM_normalize_to_bases_sequenced: false,
        assay_for_matrix: 'BaseCounts',
        use_letters_as_shapes: false,
        sizeby: '',
        connectby: '',
        connection_orderby: '',
        textby: '',
        ellipseby: '',
        dotsize: 2,
        log2tran: false,
        tsne_perplx: '',
        max_neighbors: 15,
        permanova: false,
        plotcentroids: false,
        highlight_centroids: false,
        show_centroid_distances: false,
        calculate_centroid_distances_in_all_dimensions: false,
        addtit: '',
        cdict: '',
        grid: false,
        forceaspectratio: '',
        numthreads: 1,
        return_coordinates_matrix: false,
        permanova_permutations: 10000,
        include_components_variance_plot: false,
        class_to_ignore: 'N_A',
    });
    // Parameter display names
    const displayNames = {
        glomby: 'Glom By',
        algorithm: 'Algorithm',
        PCA_Components: 'PCA Components',
        distmethod: 'Distance Method',
        subsetby: 'Subset By',
        compareby: 'Compare By',
        colourby: 'Colour By',
        colorby: 'Color By',
        shapeby: 'Shape By',
        samplesToKeep: 'Samples to Keep',
        featuresToKeep: 'Features to Keep',
        samplesToHighlight: 'Samples to Highlight',
        ignoreunclassified: 'Ignore Unclassified',
        applyfilters: 'Apply Filters',
        featcutoff: 'Feature Cutoff',
        GenomeCompletenessCutoff: 'Genome Completeness Cutoff',
        discard_SDoverMean_below: 'Discard SD over Mean Below',
        asPPM: 'As PPM',
        normalization: 'Normalization',
        PPM_normalize_to_bases_sequenced: 'PPM Normalize to Bases Sequenced',
        assay_for_matrix: 'Assay for Matrix',
        use_letters_as_shapes: 'Use Letters as Shapes',
        sizeby: 'Size By',
        connectby: 'Connect By',
        connection_orderby: 'Connection Order By',
        textby: 'Text By',
        ellipseby: 'Ellipse By',
        dotsize: 'Dot Size',
        log2tran: 'Log2 Transform',
        tsne_perplx: 't-SNE Perplexity',
        max_neighbors: 'Max Neighbors',
        permanova: 'PERMANOVA',
        plotcentroids: 'Plot Centroids',
        highlight_centroids: 'Highlight Centroids',
        show_centroid_distances: 'Show Centroid Distances',
        calculate_centroid_distances_in_all_dimensions: 'Calculate Centroid Distances in All Dimensions',
        addtit: 'Add Title',
        cdict: 'Color Dictionary',
        grid: 'Grid',
        forceaspectratio: 'Force Aspect Ratio',
        numthreads: 'Number of Threads',
        return_coordinates_matrix: 'Return Coordinates Matrix',
        permanova_permutations: 'PERMANOVA Permutations',
        include_components_variance_plot: 'Include Components Variance Plot',
        class_to_ignore: 'Class to Ignore' 
    };

    const [ordinationData, setOrdinationData] = useState(null);
    const [objects, setObjects] = useState([]);
    const [filePath, setFilePath] = useState('');
    const [selectedObj, setSelectedObj] = useState('');

    const handleChange = (e) => {
        const { name, value, type, checked } = e.target;
        setParameters({
            ...parameters,
            [name]: type === 'checkbox' ? checked: value
        });
    };

    // Function to handle the submission of user-defined parameters
    const handleSubmit = async (e) => {
        e.preventDefault();
        console.log('Selected ExpObj:', selectedObj);
        try {
            // Combine the parameters with the selected objects and file path
            const params = {
                filePath,
                ExpObj: selectedObj,
                ...parameters
            };

            // Call IPC method to run ordination script
            const result = await window.electron.runOrdinationScript(params)

            // Update the ordination data with the result
            setOrdinationData(result);
        } catch (error) {
            console.error("Error generating ordination plot:", error);
        }
    };

    // Function that allows a file to be opened
    const handleFileUpload = async () => {
        try {
            const filePath = await window.electron.invoke('open-file-dialog');
            if (filePath) {
                setFilePath(filePath);
                const result = await window.electron.invoke('load-rdata-file', filePath);
                setObjects(result);

                // Set default selection to the first object, otherwise it will not default to any
                if (result.length > 0) {
                    setSelectedObj(result[0]);
                }
            }
        } catch (error) {
            console.error('Error opening file dialog:', error);
        }
    };

    // Function to handle selecting an object after RData file upload
    const handleObjSelect = (e) => {
        const value = e.target.value;
        setSelectedObj(value);
        console.log('Selected object state:', value);
    };

    // Function to open the generated PDF file
    const handleDownloadClick = () => {
        // Send IPC event to open the ordination plot file
        window.electron.send('open-ordination-location');
    };

    // Set default selection if objects array changes
    // Ensure selectedObj is set to the first element of objects array if not already set (this initializes the selected object when the objects array is first populated)
    useEffect(() => {
        if (objects.length > 0 && !selectedObj) {
            setSelectedObj(objects[0]);
        }
    }, [objects, selectedObj]);

    return (
        <div>
            <div style={{ position: 'absolute', top: '10px', right: '10px' }}>
                <button onClick={handleNavigateTo('home')}>
                    Go Back to Home Page
                </button>
            </div>

            <h1>Generate Ordination Plot</h1>
        <div>
            {/* File upload for RData File */}
            <h3>Upload RData File for Ordination Plot</h3>
            <Button
                component="label"
                variant="contained"
                startIcon={<CloudUploadIcon />}
                onClick={handleFileUpload}
            >
                Upload RData File
            </Button>
        </div>

        {/* Dropdown for selecting Summarized Experiment Object */}
        {objects.length > 0 && (
            <div>
                <h3>Select Summarized Experiment Object</h3>
                <select onChange={handleObjSelect} value={selectedObj}>
                    {objects.map((obj, index) => (
                        <option key={index} value={obj}>
                            {obj}
                        </option>
                    ))}
                </select>
            </div>
        )}

        {/* Form for Ordination Parameters */}
        <form onSubmit={handleSubmit}>
            {Object.keys(parameters).map((key) => (
                <div key={key}>
                    <label htmlFor={key}>{displayNames[key] || key}:</label>
                    {key === 'applyfilters' ? (
                        <select
                        id={key}
                        name={key}
                        value={parameters[key]}
                        onChange={handleChange}
                        >
                            <option value='none'>none</option>
                            <option value='light'>light</option>
                            <option value='moderate'>moderate</option>
                            <option value='stringent'>stringent</option>
                        </select>
                    ) : typeof parameters[key] === 'boolean' ? (
                        <input
                        type='checkbox'
                        id={key}
                        name={key}
                        checked={parameters[key]}
                        onChange={handleChange}
                        />
                    ) : (
                        <input
                        type="text"
                        id={key}
                        name={key}
                        value={parameters[key]}
                        onChange={handleChange}
                        />
                    )}
                </div>
            ))}
            <Button type='submit' variant='outlined' color='primary'>
                Generate Ordination Plot PDF
            </Button>
        </form>

        {ordinationData && (
            <div id="ordination=container">
                <h2>Ordination Plot</h2>
                <Button
                variant='contained'
                color='primary'
                onClick={handleDownloadClick}
                >
                    Open Ordination Plot PDF
                </Button>
            </div>
        )}
        </div>
    );
};

export default Ordination;
