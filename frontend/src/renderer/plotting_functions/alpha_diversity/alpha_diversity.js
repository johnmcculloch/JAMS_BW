// alpha_diversity.js

document.getElementById('home-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to home page
    window.electron.send('navigate-to', 'renderer/home/home.html');
});

// Listener to open the alpha_diversity pdf when button is clicked
document.getElementById('openAlphaDiversityLocation').addEventListener('click', () => {
    window.electron.send('open-alphadiversity-location')
})
