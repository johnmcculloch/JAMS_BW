const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');

module.exports = {
  mode: 'development',
  entry: './src/index.js', // Entry point for the application
  output: {
    path: path.resolve(__dirname, 'dist'), // Output directory
    filename: 'bundle.js', // Output bundle file
  },
  module: {
    rules: [
      {
        test: /\.css$/, // Rule for CSS files
        use: ['style-loader', 'css-loader'], // Use style-loader and css-loader
        include: path.resolve(__dirname, 'src/renderer/shared/styles'), // Ensure CSS is picked from the right directory
      },
      {
        test: /\.jsx?$/, // Rule for JavaScript and JSX files
        exclude: /node_modules/,
        use: 'babel-loader', // Use Babel for JavaScript files
      },
    ],
  },
  plugins: [
    new HtmlWebpackPlugin({
      template: './src/renderer/home/home.html', // Source HTML file
      filename: 'home.html', // Output HTML file in the dist directory
    }),
  ],
  resolve: {
    extensions: ['.js', '.jsx'], // Resolve these file extensions
  },
};
