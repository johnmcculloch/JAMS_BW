import React, { useState } from 'react'

const Heatmap = () => {
    const [parameters, setParameters] = useState({
        glomby: '',
        hmtype: 'exploratory',
        samplesToKeep: '',
        featuresToKeep: '',
        applyfilters: 'none',
        featcutoff: '',
        GenomeCompletenessCutoff: '',
        PPM_normalize_to_bases_sequenced: false,
        subsetby: '',
        compareby: '',
        invertbinaryorder: false,
        showonlypbelow: '',
        adj_pval_for_threshold: false,
        ntop: '',
        colcategories: '',
        splitcolsby: '',
        ordercolsby: '',
        textby: '',
        cluster_samples_per_heatmap: false,
        cluster_features_per_heatmap: false,
        label_samples: false,
        cluster_rows: false,
        max_rows_in_heatmap: '',
        no_underscores: false,
        showGram: false,
        show_GenomeCompleteness: false,
        addtit: '',
        hmasPA: false,
        threshPA: 0,
        cluster_column_slices: true,
        column_split_group_order: '',
        row_order: '',
        discard_SDoverMean_below: '',
        maxl2fc: '',
        minl2fc: '',
        fun_for_l2fc: 'geom_mean',
        showpval: true,
        showroundedpval: true,
        showl2fc: true,
        assay_for_matrix: 'BaseCounts',
        normalization: 'relabund',
        asPPM: true,
        scaled: false,
        cdict: '',
        maxnumheatmaps: '',
        numthreads: 1,
        statsonlog: false,
        ignoreunclassified: true,
        returnstats: false,
        class_to_ignore: 'N_A'
    });
    const [heatmapData, setHeatmapData] = useState(null);

    const handleChange = (e) => {
        const { name, value, type, checked } = e.target;
        setParameters({
            ...parameters,
            [name]: type === 'checkbox' ? checked : value
        });
    };
    const handleSubmit = (e) => {
        e.preventDefault();


        // add logic to generate heatmap based on parameters here

    };
    return (
        <div>
            <h1>Generate Heatmap</h1>
            <form onSubmit={handleSubmit}>
                {Object.keys(parameters).map((key) => (
                    <div key={key}>
                        <label htmlFor={key}>{key}:</label>
                        {typeof parameters[key] === 'boolean' ? (
                            <input
                            type="checkbox"
                            id={key}
                            name={key}
                            checked={parameters[key]}
                            onChange={handleChange}
                            />
                        ): (
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
                <button type="submit">Generate Heatmap</button>
            </form>
            {heatmapData && (
                <div id="heatmap=container">
                    <h2>Heatmap</h2>
                    <p>{heatmapData}</p>
                    {/* Placeholder for actual heatmap visualization */}
                </div>
            )}
        </div>
    );
};

export default Heatmap;
