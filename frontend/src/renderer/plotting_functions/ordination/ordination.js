// ordination.js

document.getElementById('home-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/home/home.html');
});
