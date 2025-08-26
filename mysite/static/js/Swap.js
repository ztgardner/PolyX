

// Define 3D data
const example_data = true;
let x = [];
let y = [];
let z = [];
let color = [];
let x_vec_start = [];
let y_vec_start = [];
let z_vec_start = [];
let x_vec_end = [];
let y_vec_end = [];
let points = []
let vectors = [];
let jsonFilename = [];
let jsonObject = [];    //This is the JSON object that will contain all itp and cord info
let dropdownOptions = ["Upload Your Files", "For More", "Options"]; //Change to ITP input
let dataForPlot = []




//Loading in data to PLotly and Plotting
if (example_data) {
    x = [-10, -10, -11, -11, -11, -11, -12, -12, -12, -12, -13, -13, -13, -13, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -16, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -18, -19, -19, -19, -19, -19, -19, -20, -20, -20, -20, -20, -20, -21, -21, -21, -21, -21, -22, -22, -22, -22, -22, -23, -23, -23, -23, -23, -23, -24, -24, -24, -24, -24, -24, -24, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -25, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -26, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -27, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -30, -30, -30, -30, -30, -30, -30, -33, -33, -33, -33, -33, -33, -33, -33, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -35, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -36, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -37, -38, -38, -38, -38, -38, -38, -38, -38, -38, -38, -39, -39, -39, -39, -39, -39, -39, -40, -40, -40, -40, -40, -41, -41, -41, -41, -42, -42, -42, -42, -43, -43, -43, -43, -43, -44, -44, -44, -44, -44, -44, -45, -45, -45, -45, -45, -45, -45, -45, -45, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -46, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -48, -48, -48, -48, -48, -48, -48, -48, -48, -48, -48, -48, -48, -48, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -50, -50, -50, -50, -50, -50, -50, -53, -53, -53, -53, -54, -54, -54, -54, -55, -55, -55, -55, -55, -55, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -60, -60, -61, -61, -62, -63, -63, -64, -64, -64, -64, -64, -64, -65, -65, -65, -65, -65, -65, -65, -65, -66, -66, -66, -66, -66, -66, -66, -66, -66, -66, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -70, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -73, -73, -73, -73, -73, -73, -73, -73, -73, -73, -73, -73, -74, -74, -74, -74, -74, -74, -74, -75, -75, -75, -75, -75, -75, -75, -75, -76, -76, -76, -76, -76, -76, -76, -76, -77, -77, -77, -77, -77, -77, -77, -77, -77, -77, -78, -78, -78, -78, -78, -78, -78, -78, -78, -79, -79, -79, -79, -79, -79, -79, -80, -80, -80, -80, -81, -81, -82, -82, -83, -83, -83, -84, -84, -84, -84, -85, -85, -85, -85, -86, -86, -86, -86, -86, -87, -87, -87, -87, -87, -87, -87, -88, -88, -88, -88, -88, -88, -88, -88, -88, -89, -89, -89, -89, -89, -89, -89, -89, -89, -89, -89, -90, -90, -90, -90, -90, -90, -90, -90, -90, -90, -90, -90, -90, -90, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -91, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -93, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -94, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -95, -96, -96, -96, -96, -96, -96, -96, -96, -96, -96, -96, -97, -97, -97, -97, -97, -97, -97, -97, -97, -98, -98, -98, -98, -98, -98, -98, -98, -98, -98, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -102, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -104, -104, -104, -104, -104, -104, -104, -104, -104, -104, -104, -104, -104, -104, -105, -105, -105, -105, -105, -105, -105, -105, -105, -105, -105, -105, -106, -106, -106, -106, -106, -106, -106, -106, -106, -107, -107, -107, -107, -107, -107, -107, -108, -108, -108, -108, -108, -109, -109, -109, -109, -110, -110, -110, -110, -111, -111]
    y = [19, 45, 19, 20, 44, 45, 19, 20, 44, 45, 19, 20, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 32, 33, 44, 45, 19, 20, 32, 33, 44, 45, 19, 20, 32, 33, 45, 19, 20, 32, 33, 45, 19, 20, 21, 31, 32, 33, 19, 20, 21, 22, 31, 32, 33, 19, 20, 21, 22, 23, 24, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 23, 24, 25, 26, 27, 28, 29, 33, 34, 35, 36, 37, 38, 39, 40, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 28, 29, 30, 31, 40, 41, 42, 43, 44, 45, 28, 29, 30, 42, 43, 44, 45, 28, 29, 43, 44, 45, 28, 29, 44, 45, 28, 29, 44, 45, 28, 29, 30, 44, 45, 28, 29, 30, 43, 44, 45, 28, 29, 30, 31, 32, 42, 43, 44, 45, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 33, 34, 35, 36, 37, 38, 39, 19, 20, 44, 45, 19, 20, 44, 45, 19, 20, 21, 43, 44, 45, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 44, 45, 44, 45, 45, 28, 29, 28, 29, 51, 52, 53, 54, 28, 29, 30, 31, 51, 52, 53, 54, 28, 29, 30, 31, 32, 33, 51, 52, 53, 54, 28, 29, 30, 31, 32, 33, 34, 35, 51, 52, 53, 54, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 51, 52, 53, 54, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 51, 52, 53, 54, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 50, 51, 52, 53, 28, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 47, 48, 49, 50, 51, 52, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 40, 41, 42, 43, 44, 45, 46, 37, 38, 39, 40, 41, 42, 43, 44, 28, 35, 36, 37, 38, 39, 40, 41, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39, 28, 29, 30, 31, 32, 33, 34, 35, 36, 28, 29, 30, 31, 32, 33, 34, 28, 29, 30, 31, 28, 29, 28, 29, 19, 44, 45, 19, 20, 44, 45, 19, 20, 44, 45, 19, 20, 43, 44, 45, 19, 20, 21, 42, 43, 44, 45, 19, 20, 21, 22, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 37, 38, 39, 40, 41, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 36, 37, 38, 39, 40, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 35, 36, 37, 38, 39, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 35, 36, 37, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 19, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 27, 28, 29, 30, 31, 32, 33, 34, 35, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 19, 20, 25, 26, 27, 28, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 19, 20, 23, 24, 25, 26, 27, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 22, 23, 24, 25, 26, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 24, 38, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 23, 39, 40, 41, 42, 43, 44, 45, 19, 20, 21, 22, 41, 42, 43, 44, 45, 19, 20, 21, 42, 43, 44, 45, 19, 20, 43, 44, 45, 19, 20, 44, 45, 19, 20, 44, 45, 44, 45]
    z = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
    color = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054]
    color = color.reverse()
}
loadCords_Vectors()
//layout
const layout = {
    title: 'PolyX',

    margin: {
        l: 0,
        r: 0,
        b: 0,
        t: 40
    },
    paper_bgcolor: 'transparent',
    plot_bgcolor: 'transparent',

    scene: {
        xaxis: { title: 'X Axis',
            showgrid: false,
            showticklabels: false,
            zeroline: false,
            },
        yaxis: { title: 'Y Axis',
            showgrid: false,
            showticklabels: false,
            zeroline: false,
            },
        zaxis: { title: 'Z Axis',
            showgrid: false,
            showticklabels: false,
            zeroline: false,
            }
    }
};


data = [...points,...vectors];
Plotly.newPlot('plot-container', data, layout);

//Sets up the Dropdown Options (Dependent on ITP)
const dropdown = document.querySelector(".dropdown"); //dropdown var selection
setOptions();

// Handle dropdown selection
dropdown.addEventListener("change", () => {



    vec = mkVectors(jsonObject[dropdown.value]["atoms"]) //in index form
    vectToCords(vec)
    loadCords_Vectors() //Updates the new vectors

    dataForPlot = [...points, ...vectors];
    //console.log(dataForPlot)


    callClickListener()
    Plotly.newPlot('plot-container', dataForPlot, layout);
});

//Sets Drop down options to whatever is in dropdownOptions var
function setOptions() {
    dropdown.innerHTML = ""; // Clear existing options
    dropdownOptions.forEach(optionText => {
        const option = document.createElement("option");
        option.value = optionText;
        option.textContent = optionText;
        dropdown.appendChild(option);
    })
}

//Extend Button
function clickExtend(){
    document.getElementById('extend-button').addEventListener('click', function() {
      // Capture values from input fields
      const dihedral1 = document.getElementById("dihedral1").value
          .trim()
          .split(/\s+/)
          .map(Number);
      const dihedral2 = document.getElementById("dihedral2").value
          .trim()
          .split(/\s+/)
          .map(Number);
      const propagation = document.getElementById("propagation").value
          .trim()
          .split(/\s+/)
          .map(Number);
      const monNum = parseInt(document.getElementById("Mon_num").value);

      const data = {
          dihedral1: dihedral1,
          dihedral2: dihedral2,
          propagation: propagation,
          Mon_num: monNum
      };

      // Send data to Flask
      fetch('/Extend/extend_action', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify(data)
      })
      .then(response => {
        if (!response.ok) {
          throw new Error('Network response was not ok');
        }
        return response.json();
      })
      .then(data => {
        console.log("✅", data.message);

        if (data.files && Array.isArray(data.files)) {
          data.files.forEach(filename => {
            const link = document.createElement("a");
            link.href = `/Extend/upload/${filename}`;
            link.download = filename;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
          });
        } else {
          alert("Extension complete, but no files returned.");
        }
      })
      .catch((error) => {
        console.error('❌ Error:', error);
        alert('Extension failed. Check console for details.');
      });
    });
  }
//ITP
function click_itp(){
    document.getElementById('itp_input').click();
    updateStatusITP("Load")
}

function loadITP() {
    const file = document.getElementById('itp_input').files[0];  // Corrected file input ID
    if (file) {
        var formData = new FormData();
        formData.append('file', file);

        fetch('upload', {
            method: 'POST',
            body: formData
        })
        .then(response => response.text())
        .then(data => {


            if (data === "Files uploaded successfully" || data === "Please Upload Other File") {
                updateStatusITP(true);
            } else {
                updateStatusITP(false)
            }

        })
        .catch(error => {
            console.error('ITP Error:', error);

            updateStatusITP(false)
        });
    }
    loadJSON(); //Trys to load JSON on fileupload
}
//GRO
function click_gro(){
    document.getElementById('gro_input').click();
    updateStatusGRO("Load")
}

function loadGRO() {
    const file = document.getElementById('gro_input').files[0];  // Corrected file input ID
    if (file) {
        var formData = new FormData();
        formData.append('file', file);

        fetch('upload', {
            method: 'POST',
            body: formData
        })
        .then(response => response.text())
        .then(data => {
            //console.log(data);

            if (data === "Files uploaded successfully" || data === "Please Upload Other File") {

                updateStatusGRO(true);
            } else {

                updateStatusGRO(false)
            }
        })
        .catch(error => {
            console.error('GRO Error:', error);

            updateStatusGRO(false)
        });
    }
    loadJSON(); //Trys to load JSON on fileupload
}

function loadJSON() {
    //Load JSON into the obejct jsonObject
    //also triggers onJsonUpload()
    const file_itp = document.getElementById('itp_input').files[0]; // Get the .itp file
    const file_gro = document.getElementById('gro_input').files[0]; // Get the .gro file

    if (file_itp && file_gro) {
        // Construct the URL for the JSON file
        jsonFilename = 'upload/' + file_itp.name.replace(/\.itp$/, '.json');
        console.log('Fetching:', jsonFilename);

        //Backend needs a couple seconds sometimes to generate JSON
        function fetchWithRetry(url, retries = 5, delay = 1000) {
            return new Promise((resolve, reject) => {
                function attemptFetch(remainingRetries) {
                    fetch(url)
                        .then(response => {
                            if (!response.ok) {
                                throw new Error(`Failed to fetch ${url}: ${response.status} ${response.statusText}`);
                            }
                            return response.json();
                        })
                        .then(data => resolve(data))
                        .catch(error => {
                            if (remainingRetries > 0) {
                                console.warn(`Retrying... (${remainingRetries} retries left)`);
                                setTimeout(() => attemptFetch(remainingRetries - 1), delay);
                            } else {
                                reject(error);
                            }
                        });
                }
                attemptFetch(retries);
            });
        }

        // Use the fetchWithRetry function
        fetchWithRetry(jsonFilename)
            .then(data => {
                // Process the JSON data
                jsonObject = data;
                console.log('Fetched JSON data:', jsonObject);

                //Call Function to handle response on upload
                onFullUpload();

            })
            .catch(error => {
                console.error('Error fetching JSON file:', error.message);
            });
    } else {
        console.error('Please upload both .itp and .gro files.');
    }
}

//MISC Functions

function callClickListener() {
    // This needs to be called before drawing the plot to set up click functionality
    Plotly.newPlot('plot-container', dataForPlot, layout).then(() => {

        const plotContainer = document.getElementById('plot-container');
        let activeTextbox = null; // Track the currently active textbox

        // Get the textboxes and add focus listeners to track which is active
        const dihedral1Input = document.getElementById("dihedral1");
        const dihedral2Input = document.getElementById("dihedral2");
        const propagationInput = document.getElementById("propagation");

        [dihedral1Input, dihedral2Input, propagationInput].forEach((textbox) => {
            textbox.addEventListener("focus", () => {
                activeTextbox = textbox; // Set active textbox when focused
                textBoxToColor()//Read all indexes in textbox, change colors, on textbox click
            });
        });

        plotContainer.on('plotly_click', (data) => {
            if (data && data.points && data.points.length > 0) {
                const pointData = data.points[0]; // Assume one point clicked
                const atomIndex = pointData.pointNumber + 1; // ITP index is Plotly index + 1
                console.log(`Clicked point: (${pointData.x}, ${pointData.y}, ${pointData.z}), Index: ${atomIndex}`);

                // Ensure the clicked point is valid
                if (pointData.pointNumber >= 0 && pointData.pointNumber < x.length) {

                    if (activeTextbox) {
                        const currentValue = activeTextbox.value.trim();

                        // Add the index to the textbox if not already present and max 4 indices
                        if (!currentValue.includes(atomIndex) && currentValue.split(" ").length < 4) {
                            activeTextbox.value = currentValue
                                ? currentValue + " " + atomIndex
                                : atomIndex; // Append the atom index

                                textBoxToColor() //Makes sure that the old colors are around
                        } else {

                            //console.warn("Index already added or textbox is full!");
                        }
                    } else {
                        //console.warn("Please select a textbox to populate.");
                    }
                }

            } else {
                console.warn("No valid points found in click event data.");
            }
        });
    });
}

//Need a  function that handles manual textbox coloring input
function textBoxToColor() {
    const dihedral1Input = document.getElementById('dihedral1').value.split(" ");
    const dihedral2Input = document.getElementById('dihedral2').value.split(" ");
    const propagationInput = document.getElementById('propagation').value.split(" ");

    const plotData = Plotly._fullData;

    const newAtomColorIndex = [...color]
    const size = new Array(newAtomColorIndex.length).fill(10) //Creates an array to set the sizes of the atoms

    for (let index of dihedral1Input) {
        newAtomColorIndex[index - 1] = 'pink'
        size[index - 1] = 20
    }
    for (let index of dihedral2Input) {
        newAtomColorIndex[index - 1] = 'purple'
        size[index - 1] = 20
    }
    for (let index of propagationInput) {
        newAtomColorIndex[index - 1] = 'blue'
        size[index - 1] = 20
    }
    console.log(newAtomColorIndex)
    Plotly.update('plot-container', {
        'marker.color': [newAtomColorIndex], // Highlight color
        'marker.size': [size],
    });
}

//Deprecated Function Remove soon 1/28/2025
function colorAtomsFromInput(atomIndex,activeTextbox,color){
    atomIndex = atomIndex.split(" ")

    const newAtomColorIndex = [...color]
    if (activeTextbox === "dihedral1") {newColor = 'pink'}
        else if (activeTextbox === "dihedral2") {newColor = 'purple'}
            else if (activeTextbox === "propagation") {newColor = 'blue'}

    const size = new Array(newAtomColorIndex.length).fill(15) //Creates an array to set the sizes of the atoms

    for (let index of atomIndex) {
        newAtomColorIndex[index - 1] = newColor
        size[index - 1] = 20
    }


    Plotly.update('plot-container', {
        'marker.color': [newAtomColorIndex], // Highlight color
        'marker.size': [size],
    }
    );

}


function atomInfoFromIndex(index){
    //Finds the highlighted section
    try {
        section = dropdown.value;
    } catch {
        console.log('No DropDown Selction to Return Vectors')
        return
    }

    //Builds string to send to textbox

    const atoms = jsonObject[section]["atoms"];
    const relaventInfo = atoms.filter(list => list.includes(String(index +1))) ;//Add one to get ptoly index to itp
    let string  = atoms.filter(list => list.includes(String(index +1)));


    string.unshift(jsonObject[section]["Comments_top"]);
    string.push(jsonObject[section]["Comments_bottom"]);
    string = string.join("\n");




    //Builds Vectors for Relavent Vectors
    vec = mkVectors(relaventInfo)
    vectToCords(vec)
    loadCords_Vectors() //Updates the new vectors

}

//Takes a list of indexes and based o what section is populated in the drop down will return an array of vectors (still as indexes)
function mkVectors(index) {
    try {
        section = dropdown.value;
    } catch {
        section = "bond";
    }
    let pairs = [];

    if (section.includes("atoms")) {return []}


    if (section.includes("virtual site")) {
        for (let entry of index) {
            // Ensure we are iterating through the actual elements in the array
            for (let i = 0; i < entry.length - 1; i++) {
                pairs.push([entry[0], entry[i + 1]]);
            }
        }
        return pairs;
    }



    // entry is whatever itp line is, i.e., angle is [1,2,3], bond is [1,2], etc
    for (let entry of index) {


        // Ensure we are iterating through the actual elements in the array
        for (let i = 0; i < entry.length - 1; i++) {
            if (entry[i].trim() === "") {continue} //skip empty string
            pairs.push([entry[i], entry[i + 1]]);
        }
    }

    //Removes Dupes
    //pairs = Array.from(new Set(pairs.map(a => JSON.stringify(a)))).map(e => JSON.parse(e));
    return pairs;
}

//passes a list of vectors (in index form ideally from mkVec()) and converts it to a list of cordinates
//loads into the vec_start/end varaibles
function vectToCords(vect_list){

    x_vec_start = [];
    y_vec_start = [];
    z_vec_start = [];
    x_vec_end = [];
    y_vec_end = [];
    z_vec_end = [];

    for (let entry of vect_list) {
        entry = entry.map(num => parseInt(num, 10))

        //console.log(entry)


        x_vec_start.push(jsonObject.coordinates[entry[0] - 1][0]);
        y_vec_start.push(jsonObject.coordinates[entry[0] - 1][1]);
        z_vec_start.push(jsonObject.coordinates[entry[0] - 1][2]);

        x_vec_end.push(jsonObject.coordinates[entry[1] - 1][0]);
        y_vec_end.push(jsonObject.coordinates[entry[1] - 1][1]);
        z_vec_end.push(jsonObject.coordinates[entry[1] - 1][2]);


        //console.log(x_vec_start,y_vec_start,z_vec_start)
        //console.log(x_vec_end,y_vec_end,z_vec_end)
    }

}


function loadCords_Vectors() {

    points = [
        {
            x: x,
            y: y,
            z: z,
            mode: 'markers',
            marker: {
                size: 5,
                color: color,
                colorscale: 'Viridis',
                opacity: 0.8
            },
            type: 'scatter3d',

        }
        ];
    vectors = [
        {


            x: x_vec_start.flatMap((start, i) => [start, x_vec_end[i], null]), // Start, End, Null
            y: y_vec_start.flatMap((start, i) => [start, y_vec_end[i], null]),
            z: z_vec_start.flatMap((start, i) => [start, z_vec_end[i], null]),
         /*
            u: x_vec_end.map((x, i) => x - x_vec_start[i]),
            v: y_vec_end.map((y, i) => y - y_vec_start[i]),
            w: z_vec_end.map((z, i) => z - z_vec_start[i]),
         */
            mode: 'lines',
            marker: {
                size: 20,
                color: 'red',
                colorscale: 'Viridis',
                opacity: 0.8
            },
            type: 'scatter3d',
            hoverinfo: "none"
        }
    ];
}
//everyone needs fun colors//everyone needs fun colors
function elementsToColors(elements) {
    const elementColors = {
        'H': '#808080', // Hydrogen - Grey
        'C': '#000000', // Carbon - Black
        'N': '#0000FF', // Nitrogen - Blue
        'O': '#FF0000', // Oxygen - Red
        'P': '#FFA500', // Phosphorus - Orange
        'S': '#FFFF00', // Sulfur - Yellow
        'F': '#00FF00', // Fluorine - Green
        'Cl': '#00FF00', // Chlorine - Green
        'Br': '#8B4513', // Bromine - Brown
        'I': '#800080', // Iodine - Purple

    };
    return elements.map(element => elementColors[element[0].charAt(0).toUpperCase()] || '#808080'); // Default color for unknown elements (Grey)
}
//Helper function to extract the jsonObject coridnates into x,y,z cords plotly uses
function jsonToPoints(json_data) {
    //console.log(json_data);
    //console.log(json_data.coordinates);
    const coordinates = json_data.coordinates;

    x = coordinates.map(coord => coord[0]); y = coordinates.map(coord => coord[1]); z = coordinates.map(coord => coord[2]);

    points = [
    {
        x: x,
        y: y,
        z: z,
        mode: 'markers',
        marker: {
            size: 5,
            color: color,
            colorscale: 'Viridis',
            opacity: 0.8
        },
        type: 'scatter3d'
    }
    ];

}

//Handles Upload Checkmarks
function updateStatusITP(success) {
    const statusElement = document.getElementById('upload-status-itp');

    if (success === true) {
        statusElement.textContent = '✔'; // Checkmark character
        statusElement.style.color = 'green';
    } else if (success === false) {
        statusElement.textContent = '✘'; // Xmark
        statusElement.style.color = 'red';
    } else {
        statusElement.textContent = '⏳'; // Hourglass character
        statusElement.style.color = '#0033A0'; // Optional: blue for loading
    }
}

function updateStatusGRO(success) {
    const statusElement = document.getElementById('upload-status-gro');

    if (success === true) {
        statusElement.textContent = '✔'; // Checkmark character
        statusElement.style.color = 'green';
    } else if (success === false) {
        statusElement.textContent = '✘'; // Xmark
        statusElement.style.color = 'red';
    } else {
        statusElement.textContent = '⏳'; // Hourglass character
        statusElement.style.color = '#0033A0'; // Optional: blue for loading
    }
}

//inputs a index or list of indexes, replaces index with corindates of that atom




let itpSection = []
function onFullUpload() {



    //updates drop down to itp sections
    itpSection =  Object.keys(jsonObject)
    itpSection = itpSection.filter(name => !["moleculetype", "coordinates"].includes(name));
    dropdownOptions = itpSection
    setOptions()


    //console.log(jsonObject.atoms.atom_name[0]);
    color = elementsToColors(jsonObject.atoms.atom_name);


    //Populates plot with ITP data

    layout.title = jsonFilename.split('/').pop().split('.')[0] //Changes plot title to json name


    jsonToPoints(jsonObject)
    // Update the Plot



    //Listening to Clicks
    callClickListener();


    dataForPlot = [...points, ...vectors];
    Plotly.newPlot('plot-container', dataForPlot, layout);

    //Changes aspect ratio to stop weird stretching
    const xRange = Math.max(...x) - Math.min(...x);
    const yRange = Math.max(...y) - Math.min(...y);
    const zRange = Math.max(...z) - Math.min(...z);
    Plotly.relayout('plot-container', {
    'scene.aspectratio': { x: xRange, y: yRange, z: zRange }});



}


//BUILDING OFF GUI, POLYMER EXTENSION STUFF


document.addEventListener("DOMContentLoaded", () => {
    let activeTextbox = null; // Variable to track the active textbox

    // Get the textboxes
    const dihedral1Input = document.getElementById("dihedral1");
    const dihedral2Input = document.getElementById("dihedral2");
    const propagationInput = document.getElementById("propagation");

    // Add event listeners to track which textbox is active
    [dihedral1Input, dihedral2Input, propagationInput].forEach((textbox) => {
      textbox.addEventListener("focus", () => {
        activeTextbox = textbox; // Set the active textbox when focused
      });
    });

    // Example: Assuming each atom is a clickable button with a data attribute for the index
    const atomElements = document.querySelectorAll("[data-atom-index]");

    atomElements.forEach((atom) => {
      atom.addEventListener("click", () => {
        if (activeTextbox) {
          const atomIndex = atom.getAttribute("data-atom-index");
          const currentValue = activeTextbox.value.trim();

          // Ensure no duplicate entries and max 4 indices
          if (!currentValue.includes(atomIndex) && currentValue.split(" ").length < 4) {
            activeTextbox.value = currentValue
              ? currentValue + " " + atomIndex
              : atomIndex; // Append the atom index
          } else {
            alert("Index already added or textbox is full!");
          }
        } else {
          alert("Please select a textbox first!");
        }
      });
    });
  });



/* HELP ICON */
// Get the help icon and overlay
const helpIcon = document.getElementById('help-icon');
const overlay = document.getElementById('overlay');
const closeBtn = document.getElementById('close-btn');

// Show overlay when the help icon is clicked
helpIcon.addEventListener('click', function() {
    overlay.style.visibility = 'visible';
});

// Close the overlay when the close button is clicked
closeBtn.addEventListener('click', function() {
    overlay.style.visibility = 'hidden';
});

document.getElementById("help-icon").addEventListener("click", () => {
    const overlay = document.getElementById("overlay");
    overlay.classList.add("visible"); // Show overlay
});

document.getElementById("close-btn").addEventListener("click", () => {
    const overlay = document.getElementById("overlay");
    overlay.classList.remove("visible"); // Hide overlay
});