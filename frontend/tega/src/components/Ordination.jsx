import React, { useEffect, useState } from 'react';


const Ordination = () => {
    const [parameters, setParameters] = useState({
        glomby: '',
        algorithm: '',
        PCA_Components: '',
        distmehod: '',
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
        dotsize: '',
        log2tran: false,
        tsne_perplx: '',
        max_neighbors: '',
        permanova: false,
        plotcentroids: false,
        highlight_centroids: false,
        show_centroid_distances: false,
        calculate_centroid_distances_in_all_dimensions: false,
        addtit: '',
        cdict: '',
        grid: false,
        forceaspectratio: '',
        numthreads: '',
        return_coordinates_matrix: false,
        permanova_permutations: '',
        include_components_variance_plot: false,
        class_to_ignore: 'N_A',
    });

    const [objects, setObjects] = useState([]);
    const [filePath, setFilePath] = useState('');
    const [selectedObj, setSelectedObj] = useState('');

    const handleChange = (e) => {
        const { name, value, type, checked } = e.target;
        setParameters({
            ...parameters,
            [name]: type === 'checkbox' ? checked : value
        });
    };

    const handleClick = () => {
        window.electron.send('navigate-to', 'home');
    };

    return (
        <div>
            <div style={{ position: 'absolute', top: '10px', right: '10px' }}>
                <button onClick={handleClick}>
                    Go Back to Home Page
                </button>
            </div>
        </div>
    );

};

export default Ordination;
