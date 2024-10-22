module.exports = [
  // Add support for native node modules
  {
    // Support for native modules
    test: /native_modules[/\\].+\.node$/,
    use: 'node-loader',
  },
  {
    test: /[/\\]node_modules[/\\].+\.(m?js|node)$/,
    parser: { amd: false },
    use: {
      loader: '@vercel/webpack-asset-relocator-loader',
      options: {
        outputAssetBase: 'native_modules',
      },
    },
  },
  {
    // JavaScript and JSX files
    test: /\.jsx?$/,
    exclude: /node_modules/,  // Moved outside of options
    use: {
      loader: 'babel-loader',
      options: {
        presets: ['@babel/preset-env', '@babel/preset-react'], // Make sure to include @babel/preset-env
      },
    },
  },
  {
    // CSS Loader
    test: /\.css$/,
    use: ['style-loader', 'css-loader'],
  },
  {
    // File loader for images and other assets
    test: /\.(png|jpg|gif|svg)$/,
    use: ['file-loader'],
  },
];
