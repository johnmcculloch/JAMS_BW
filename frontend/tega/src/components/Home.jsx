import React from 'react'

export default function Home() {
    const handleClick = () => {
        // send IPC message to main process to navigate to heatmap page
        window.electron.send('navigate-to', 'heatmap');
    };

    return (
        <div>
            Home<br></br>
            <button onClick={handleClick}>
                Go to Relabund Heatmap Analysis
            </button>
        </div>
    );
}
