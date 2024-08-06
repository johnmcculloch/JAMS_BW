// home.js

document.getElementById('heatmap-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/plotting_functions/heatmap/heatmap.html');
});


document.getElementById('ordination-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/plotting_functions/ordination/ordination.html');
});

document.getElementById('relabundfeatures-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/plotting_functions/relabund_features/relabund_features.html');
});

document.getElementById('alphadiversity-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/plotting_functions/alpha_diversity/alpha_diversity.html');
});
