const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const path = require('path');
const { exec } = require('child_process');
const fs = require('fs');

let mainWindow;

console.log(`app.isPackaged: ${app.isPackaged}`); // debugging
let isDev = !app.isPackaged;

console.log(`isDev: ${isDev}`);

function createWindow() {
  mainWindow = new BrowserWindow({
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

  mainWindow.on('closed', function () {
    mainWindow = null;
  });

  // Check if JAMS is installed
  checkJAMSInstallation();
}

app.on('ready', createWindow);

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});

function checkJAMSInstallation() {
  const homeDir = process.env.HOME || process.env.HOMEPATH || process.env.USERPROFILE;
  const binDir = path.join(homeDir, 'bin');
  const executables = ['JAMSalpha', 'JAMSbeta', 'JAMS16', 'JAMSbuildk2db', 'JAMSjoinlanes', 'JAMSmakeswarm', 'JAMSbankit', 'JAMSfastqprefixrenamer'];

  let allExecutablesExist = true;

  for (const exec of executables) {
    const execPath = path.join(binDir, exec);
    if (!fs.existsSync(execPath)) {
      allExecutablesExist = false;
      break;
    }
  }

  if (allExecutablesExist) {
    console.log('JAMS is installed.');
  } else {
    console.log('JAMS is not installed.');
    // JAMS is not installed, prompt the user to install it
    dialog.showMessageBox({
      type: 'question',
      buttons: ['Install', 'Cancel'],
      defaultId: 0,
      title: 'Install JAMS',
      message: 'JAMS is not installed. Would you like to install it now?',
    }).then(result => {
      if (result.response === 0) {
        // User chose to install JAMS
        installJAMS();
      }
    });
  }
}

function installJAMS() {
  const installerPath = path.join(process.resourcesPath, 'JAMSinstaller_BW');
  const homeDir = process.env.HOME || process.env.HOMEPATH || process.env.USERPROFILE;
  const binDir = path.join(homeDir, 'bin');

  // Command to create the bin directory if it doesn't exist
  const createBinDirCommand = `mkdir -p ${binDir}`;

  // Command to open a new terminal window and run the installer
  const installCommand = `osascript -e 'tell application "Terminal" to do script "cd ${path.join(process.resourcesPath, 'JAMS_BW_dev')} && ${createBinDirCommand} && ./${path.basename(installerPath)} --install"'`;

  exec(installCommand, (error, stdout, stderr) => {
    if (error) {
      console.error(`Error installing JAMS: ${error.message}`);
      return;
    }
    if (stderr) {
      console.error(`Stderr: ${stderr}`);
      return;
    }
    console.log(`JAMS installed successfully: ${stdout}`);
  });
}

// ** Main Processes ** 

// Handle Page Navigation IPC event
ipcMain.on('navigate-to', (event, page) => {
  // Send the navigation event to the renderer
  mainWindow.webContents.send('navigate-to', page);
});


// Function to handle RData session files
ipcMain.handle('open-file-dialog', async () => {
  const { canceled, filePaths } = await dialog.showOpenDialog({
    properties: ['openFile'],
    filters: [
      { name: 'RData Files', extensions: ['rdata', 'rda', 'rds'] }
    ]
  });

  if (canceled) {
    return null;
  } else {
    return filePaths[0];
  }
});


// Load user-provided R data file containing Expvec
ipcMain.handle('load-rdata-file', async (event, filePath) => {
  return new Promise((resolve, reject) => {
    const rscriptPath = '/usr/local/bin/Rscript'; // Locate R on user's system

    // Check if the file is an .rds file
    const isRdsFile = path.extname(filePath).toLowerCase() === '.rds';

    const command = `
      ${rscriptPath} -e "
      file_path <- '${filePath}'
      if (${isRdsFile ? 'TRUE' : 'FALSE'}) {
        obj <- readRDS(file_path)
      } else {
        load(file_path)
      }
      list_objs <- ls()
      list_objs <- list_objs[sapply(list_objs, function(x) is.list(get(x)))]
      if (length(list_objs) > 0) {
        first_list <- get(list_objs[1])
        obj_names <- names(first_list)
        writeLines(paste(list_objs[1], obj_names, sep='$'), stdout())
      } else {
        writeLines('', stdout())
      }
      "
    `;

    exec(command, (error, stdout, stderr) => {
      if (error) {
        console.error(`exec error: ${error}`);
        reject(`Error: ${error.message}`);
        return;
      }
      if (stderr) {
        console.error(`stderr: ${stderr}`);
        if (!stderr.toLowerCase().includes("warning")) {
          reject(`R error: ${stderr}`);
          return;
        }
      }
      const objects = stdout.split('\n').filter(name => name.trim() !== '');
      resolve(objects); // Ensure objects are returned here
    });
  });
});



// ** Heatmap Script **
ipcMain.handle('run-heatmap-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Received params:', params);

  const { filePath, ExpObj, advancedSettings, ...otherParams } = params;

  // Split ExpObj to capture the file name dynamically
  const [fileName, objName] = ExpObj.split('$');
  console.log(`File Name: ${fileName}, Object Name: ${objName}`);

  // Construct paramStr dynamically to account for anything the user inputs
  const paramStr = `ExpObj = ${ExpObj}, ` +
    Object.entries(otherParams)
      .map(([key, value]) => {
        if (value === "" || value === null || value === 'null' || value === 'NULL') {
          return `${key}=NULL`;
        } else if (typeof value === 'string' && value.startsWith('c(')) { // Check if param is a R variable
          return `${key}=${value}`;
        } else if (typeof value === 'boolean') {
          return `${key}=${value ? 'TRUE' : 'FALSE'}`;
        } else if (typeof value ==='number' || !isNaN(value)) {
          return `${key}=${Number(value)}`
        } else if (typeof value === 'string') {
          return `${key}="${value}"`;
        } else {
          return `${key}=${value}`;
        }
      })

      .join(', ');

  console.log(paramStr);

  // Send paramStr to the renderer process for debugging
  event.sender.send('param-str', paramStr);

  // Determine appropriate R command based on file extension
  const fileExtension = path.extname(filePath).toLowerCase();
  const loadCommand = fileExtension === '.rds' 
    ? `obj <- readRDS("${filePath}")` 
    : `load("${filePath}")`;

// Run plot_relabund_heatmap command with user-defined parameters
  const outputDir = path.join(app.getPath('userData'), 'assets');
  const outputFilePath = path.join(outputDir, 'heatmap.pdf');
  const rscriptPath = '/usr/local/bin/Rscript';
  const scriptPath = isDev
      ? path.join(__dirname, '..', 'R', 'plot_relabund_heatmap.R') // dev path
      : path.join(process.resourcesPath, 'JAMS_BW_dev', 'R', 'plot_relabund_heatmap.R'); // prod path

  // Create output directory and log results
  try {
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`Created directory: ${outputDir}`);
    }
    console.log(`Output directory exists: ${outputDir}`);
  } catch (err) {
    console.error(`Error creating directory: ${err}`);
  }

  // Escape spaces in file path for R
  const escapedOutputPath = outputFilePath.replace(/ /g, '\\ ');

  const script = `
   ${rscriptPath} -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      options(encoding = "UTF-8");
      ${loadCommand}
      library(JAMS); 
      source("${scriptPath}");
      pdf("${escapedOutputPath}");
      tryCatch({
        plot_relabund_heatmap(${paramStr})
      }, error = function(e) {
        print(paste("Error:", e$message))
      })
      dev.off();
    })
    })'
  `;

  console.log('Executing command:', script);

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
      resolve({ stdout: stdout, imagePath: outputFilePath });
    });
  });
});

// IPC handler for opening the heatmap file
ipcMain.on('open-heatmap-location', (event, filePath) => {
  const outputFilePath = path.join(app.getPath('userData'), 'assets', 'heatmap.pdf');
  shell.openPath(outputFilePath);
});


// ** Ordination Script **
ipcMain.handle('run-ordination-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Received params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Split ExpObj to capture the file name dynamically
  const [fileName, objName] = ExpObj.split('$');
  console.log(`File Name: ${fileName}, Object Name: ${objName}`);

  // Construct paramStr dynamically to account for anything the user inputs
  const paramStr = `ExpObj = ${ExpObj}, ` +
    Object.entries(otherParams)
      .map(([key, value]) => {
        if (value === "" || value === null || value === 'null' || value === 'NULL') {
          return `${key}=NULL`;
        } else if (typeof value === 'string' && value.startsWith('c(')) { // Check if param is a R variable
          return `${key}=${value}`;
        } else if (typeof value === 'boolean') {
          return `${key}=${value ? 'TRUE' : 'FALSE'}`;
        } else if (typeof value ==='number' || !isNaN(value)) {
          return `${key}=${Number(value)}`
        } else if (typeof value === 'string') {
          return `${key}="${value}"`;
        } else {
          return `${key}=${value}`;
        }
      })

      .join(', ');

  console.log(paramStr);

  // Send paramStr to the renderer process for debugging
  event.sender.send('param-str', paramStr);

  // Determine appropriate R command based on file extension
  const fileExtension = path.extname(filePath).toLowerCase();
  const loadCommand = fileExtension === '.rds' 
    ? `obj <- readRDS("${filePath}")` 
    : `load("${filePath}")`;

// Run plot_Ordination command with user-defined parameters
  const outputDir = path.join(app.getPath('userData'), 'assets');
  const outputFilePath = path.join(outputDir, 'ordination.pdf');
  const rscriptPath = '/usr/local/bin/Rscript';
  const scriptPath = isDev
      ? path.join(__dirname, '..', 'R', 'plot_Ordination.R') // dev path
      : path.join(process.resourcesPath, 'JAMS_BW_dev', 'R', 'plot_Ordination.R'); // prod path

  // Create output directory and log results
  try {
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`Created directory: ${outputDir}`);
    }
    console.log(`Output directory exists: ${outputDir}`);
  } catch (err) {
    console.error(`Error creating directory: ${err}`);
  }

  // Escape spaces in file path for R
  const escapedOutputPath = outputFilePath.replace(/ /g, '\\ ');

  const script = `
   ${rscriptPath} -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      options(encoding = "UTF-8");
      ${loadCommand}
      library(JAMS); 
      source("${scriptPath}");
      pdf("${escapedOutputPath}");
      tryCatch({
        print(plot_Ordination(${paramStr}))
      }, error = function(e) {
        print(paste("Error:", e$message))
      })
      dev.off();
    })
    })'
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
      resolve({ stdout: stdout, imagePath: outputFilePath });
    });
  });
});


// IPC handler for opening the ordination file
ipcMain.on('open-ordination-location', (event, filePath) => {
  const outputFilePath = path.join(app.getPath('userData'), 'assets', 'ordination.pdf');
  shell.openPath(outputFilePath);
});




// ** AlphaDiversity Script **
ipcMain.handle('run-alphaDiversity-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Received params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
  const paramStr = `ExpObj = ${ExpObj}, ` +
    Object.entries(otherParams)
      .map(([key, value]) => {
        if (value === "" || value === null || value === 'null' || value === 'NULL') {
          return `${key}=NULL`;
        } else if (typeof value === 'string' && value.startsWith('c(')) { // Check if param is a R variable
          return `${key}=${value}`;
        } else if (typeof value === 'boolean') {
          return `${key}=${value ? 'TRUE' : 'FALSE'}`;
        } else if (typeof value ==='number' || !isNaN(value)) {
          return `${key}=${Number(value)}`
        } else if (typeof value === 'string') {
          return `${key}="${value}"`;
        } else {
          return `${key}=${value}`;
        }
      })

      .join(', ');

  console.log(paramStr);

  // Send paramStr to the renderer process for debugging
  event.sender.send('param-str', paramStr);

  // Determine appropriate R command based on file extension
  const fileExtension = path.extname(filePath).toLowerCase();
  const loadCommand = fileExtension === '.rds' ? `obj <- readRDS("${filePath}")` : `load("${filePath}")`;

// Run plot_alphadiversity command with user-defined parameters
  const outputDir = path.join(app.getPath('userData'), 'assets');
  const outputFilePath = path.join(outputDir, 'alphaDiversity.pdf');
  const rscriptPath = '/usr/local/bin/Rscript';
  const scriptPath = isDev
      ? path.join(__dirname, '..', 'R', 'plot_alpha_diversity.R') // dev path
      : path.join(process.resourcesPath, 'JAMS_BW_dev', 'R', 'plot_alpha_diversity.R'); // prod path

  // Create output directory and log results
  try {
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`Created directory: ${outputDir}`);
    }
    console.log(`Output directory exists: ${outputDir}`);
  } catch (err) {
    console.error(`Error creating directory: ${err}`);
  }

  // Escape spaces in file path for R
  const escapedOutputPath = outputFilePath.replace(/ /g, '\\ ');

  const script = `
    ${rscriptPath} -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      ${loadCommand};
      library(JAMS); 
      source("${scriptPath}");
      pdf("${escapedOutputPath}", paper = "a4r");
      print(plot_alpha_diversity(${paramStr}))
      dev.off();
      })
    })'
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
      resolve({ stdout: stdout, imagePath: outputFilePath });
    });
  });
});

// IPC handler for opening the alpha diversity file
ipcMain.on('open-alphadiversity-location', (event, filePath) => {
  const outputFilePath = path.join(app.getPath('userData'), 'assets', 'alphaDiversity.pdf');
  shell.openPath(outputFilePath);
});

// ** RelabundFeatures Script **
ipcMain.handle('run-relabundFeatures-script', async (event, params) => {
  // Log the params object for debugging
  console.log('Received params:', params);

  const { filePath, ExpObj, ...otherParams } = params;

  // Construct paramStr dynamically to account for anything the user inputs
  const paramStr = `ExpObj = ${ExpObj}, ` +
    Object.entries(otherParams)
      .map(([key, value]) => {
        if (value === "" || value === null || value === 'null' || value === 'NULL') {
          return `${key}=NULL`;
        } else if (typeof value === 'string' && value.startsWith('c(')) { // Check if param is a R variable
          return `${key}=${value}`;
        } else if (typeof value === 'boolean') {
          return `${key}=${value ? 'TRUE' : 'FALSE'}`;
        } else if (typeof value ==='number' || !isNaN(value)) {
          return `${key}=${Number(value)}`
        } else if (typeof value === 'string') {
          return `${key}="${value}"`;
        } else {
          return `${key}=${value}`;
        }
      })

      .join(', ');

  console.log(paramStr);

  // Send paramStr to the renderer process for debugging
  event.sender.send('param-str', paramStr);

  // Determine appropriate R command based on file extension
  const fileExtension = path.extname(filePath).toLowerCase();
  const loadCommand = fileExtension === '.rds' ? `obj <- readRDS("${filePath}")` : `load("${filePath}")`;

// Run plot_relabundfeatures command with user-defined parameters
  const outputDir = path.join(app.getPath('userData'), 'assets');
  const outputFilePath = path.join(outputDir, 'relabundFeatures.pdf');
  const rscriptPath = '/usr/local/bin/Rscript';
  const scriptPath = isDev
      ? path.join(__dirname, '..', 'R', 'plot_relabund_features.R') // dev mode
      : path.join(process.resourcesPath, 'JAMS_BW_dev', 'R', 'plot_relabund_features.R'); // prod mode


  // Create output directory and log results
  try {
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
      console.log(`Created directory: ${outputDir}`);
    }
    console.log(`Output directory exists: ${outputDir}`);
  } catch (err) {
    console.error(`Error creating directory: ${err}`);
  }

  // Escape spaces in file path for R
  const escapedOutputPath = outputFilePath.replace(/ /g, '\\ ');


  const script = `
    ${rscriptPath} -e '
    suppressPackageStartupMessages({
    suppressWarnings({
      ${loadCommand};
      library(JAMS);  
      source("${scriptPath}");
      pdf("${escapedOutputPath}", paper = "a4r");
      print(plot_relabund_features(${paramStr}))
      dev.off();
      })
    })'
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
      resolve({ stdout: stdout, imagePath: outputFilePath });
    });
  });
});


// IPC handler for opening the relabund features file
ipcMain.on('open-RelabundFeatures-location', (event, filePath) => {
  const outputFilePath = path.join(app.getPath('userData'), 'assets', 'relabundFeatures.pdf');
  shell.openPath(outputFilePath);
});
