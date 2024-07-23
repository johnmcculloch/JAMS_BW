const { app, BrowserWindow, ipcMain } = require('electron');
const { exec } = require('child_process');
const path = require('path');

// Handle creating/removing shortcuts on Windows when installing/uninstalling.
if (require('electron-squirrel-startup')) {
  app.quit();
}

const createWindow = () => {
  // Create the browser window.
  const mainWindow = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      enableRemoteModule: false,
      nodeIntegration: false,
    },
  });

  // and load the index.html of the app.
  mainWindow.loadFile(path.join(__dirname, 'index.html'));


  // Open the DevTools.
  mainWindow.webContents.openDevTools();
};

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.whenReady().then(() => {
  createWindow();


  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});


app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

// R Heatmap script

ipcMain.handle('load-rdata-file', async (event, filePath) => {
  const script = `
    Rscript -e '
      load("${filePath}");
      objs <- names(expvec_MR50);
      writeLines(objs, stdout());
    '
  `;

  return new Promise((resolve, reject) => {
    exec(script, (error, stdout, stderr) => {
      if (error) {
        reject(`Error: ${error.message}`);
        return;
      }
      if (stderr) {
        reject(`Stderr: ${stderr}`);
        return;
      }
      const objects = stdout.split('\n').filter(name => name.trim() !== '');
      resolve(objects); // Ensure objects are returned here
    });
  });
});



ipcMain.handle('run-r-script', async (event, params) => {
  const paramStr = Object.entries(params)
    .map(([key, value]) => `${key}=${typeof value === 'string' ? `"${value}"` : value}`)
    .join(', ');

  const script = `
    Rscript -e 'library(JAMS); source("/Users/mossingtonta/Projects/JAMS_BW/R/plot_relabund_heatmap.R"); plot_relabund_heatmap(${paramStr})'
  `;

  return new Promise((resolve, reject) => {
    exec(script, (error, stdout, stderr) => {
      if (error) {
        reject(`Error: ${error.message}`);
        return;
      }
      if (stderr) {
        reject(`Stderr: ${stderr}`);
        return;
      }
      resolve(`Stdout: ${stdout}`);
    });
  });
});
