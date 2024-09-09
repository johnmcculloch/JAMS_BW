/******/ (() => { // webpackBootstrap
/******/ 	var __webpack_modules__ = ({

/***/ "electron":
/*!***************************!*\
  !*** external "electron" ***!
  \***************************/
/***/ ((module) => {

"use strict";
module.exports = require("electron");

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/compat */
/******/ 	
/******/ 	if (typeof __webpack_require__ !== 'undefined') __webpack_require__.ab = __dirname + "/native_modules/";
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
/*!************************!*\
  !*** ./src/preload.js ***!
  \************************/
var _require = __webpack_require__(/*! electron */ "electron"),
  contextBridge = _require.contextBridge,
  ipcRenderer = _require.ipcRenderer;
contextBridge.exposeInMainWorld('electron', {
  send: function send(channel, data) {
    return ipcRenderer.send(channel, data);
  },
  receive: function receive(channel, func) {
    return ipcRenderer.on(channel, function (event) {
      for (var _len = arguments.length, args = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
        args[_key - 1] = arguments[_key];
      }
      return func.apply(void 0, args);
    });
  },
  loadRDataFile: function loadRDataFile(filePath) {
    return ipcRenderer.invoke('load-rdata-file', filePath);
  },
  invoke: function invoke(channel, data) {
    return ipcRenderer.invoke(channel, data);
  }
});
module.exports = __webpack_exports__;
/******/ })()
;
//# sourceMappingURL=preload.bundle.js.map