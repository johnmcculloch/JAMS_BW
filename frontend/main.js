const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const path = require('path');

console.log(`app.isPackaged: ${app.isPackaged}`); // Add this line for debugging
let isDev = !app.isPackaged;

console.log(`isDev: ${isDev}`);

function createWindow() {
  const mainWindow = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
    },
  });

  const startUrl = isDev
    ? 'http://localhost:1234'
    : `file://${path.join(__dirname, 'dist', 'index.html')}`;
  console.log(`Loading URL: ${startUrl}`); 
  mainWindow.loadURL(startUrl);
}

app.on('ready', createWindow);

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

// ** Main Processes ** 

// Handle Page Navigation IPC event
ipcMain.on('navigate-to', (event, page) => {
  // Send the navigation event to the renderer
  mainWindow.webContents.send('navigate-to', page);
});
