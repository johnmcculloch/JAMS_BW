import React from 'react';
import Button from '@mui/material/Button';

export default function Home() {
    const handleClick = () => {
        // send IPC message to main process to navigate to heatmap page
        window.electron.send('navigate-to', 'heatmap');
    };

    return (
        <div>
            Home<br></br>
            <Button
                    variant="outlined"
                    color='primary'
                    onClick={handleClick}
                    >
                        Go to Relabund Heatmap Analysis
                    </Button>
        </div>
    );
}
