import React from 'react';
import Button from '@mui/material/Button';

export default function Home() {
    const handleNavigateToHeatmap = () => {
        // send IPC message to main process to navigate to heatmap page
        window.electron.send('navigate-to', 'heatmap');
    };

    const handleNavigateToOrdination = () => {
        // send IPC message to main process to navigate to ordination page
        window.electron.send('navigate-to', 'ordination');
    };

    const handleNavigateToAlphaDiversity = () => {
        // Send IPC message to main process to navigate to alpha diversity page
        window.electron.send('navigate-to', 'alphadiversity')
    };

    return (
        <div>
            Home<br></br>
            <Button
                    variant="outlined"
                    color='primary'
                    onClick={handleNavigateToHeatmap}
                    >
                        Go to Relabund Heatmap Analysis
                    </Button>
            
            <Button
                variant="outlined"
                color='primary'
                onClick={handleNavigateToOrdination}
                >
                    Go to Plot Ordination Analysis
                </Button>

            <Button
                variant='outlined'
                color='primary'
                onClick={handleNavigateToAlphaDiversity}
                >
                    Go to Alpha Diversity Analysis
                </Button>
        </div>
    );
}
