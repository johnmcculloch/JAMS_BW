// home.js

document.getElementById('heatmap-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/plotting_functions/heatmap/heatmap.html');
});

