import { useState, useEffect, useRef, useCallback } from 'react';
import { FaTrash, FaDownload, FaGlobe, FaUpload, FaLightbulb, FaRandom, FaBrain, FaCut, FaMagic, FaChevronUp } from 'react-icons/fa';
import './App.css';

// ============================================================================
// PERFORMANCE LIMITS - Keep large graphs from freezing the UI
// ============================================================================

// Hard caps to keep algorithms from doing too much work on huge graphs.
// These values are conservative and aimed at keeping the browser responsive.
// Aggressively cap crossing comparisons to prevent noticeable freezes even on
// mid-sized dense graphs. We only need an approximate crossing count for UX.
const MAX_CROSSING_COMPARISONS = 10_000;       // Max edge-pair checks when counting crossings
const MAX_NODES_FOR_PLANARITY = 250;          // Above this, skip NetworkX planarity checks
const MAX_EDGES_FOR_PLANARITY = 800;

// NOTE: The minimum node/edge removal search is *combinatorial* and can explode
// even for dense medium-sized graphs. To avoid the UI freezing when loading
// JSON graphs (e.g. ~15 nodes but very dense), we only run these searches on
// very small graphs.
const MAX_NODES_FOR_MIN_REMOVAL = 20;         // Above this, skip expensive min-removal searches
const MAX_EDGES_FOR_MIN_REMOVAL = 40;

// Additional guardrails specifically for uploaded graphs. Rendering and even
// basic analysis of extremely large graphs can freeze the browser, so we
// reject files that exceed these limits with a clear message instead.
// These limits are intentionally conservative to favor responsiveness.
const MAX_UPLOAD_NODES = 400;
const MAX_UPLOAD_EDGES = 1200;

// ============================================================================
// PYODIDE / NETWORKX INTEGRATION - For proper planarity testing and layout
// ============================================================================

let pyodideInstance = null;
let pyodideLoading = false;
let pyodideLoadPromise = null;

/**
 * Loads Pyodide and NetworkX for use in planarity testing and layout.
 * Uses lazy loading - only loads when first needed.
 */
async function loadPyodideAndNetworkX() {
  if (pyodideInstance) return pyodideInstance;
  if (pyodideLoadPromise) return pyodideLoadPromise;

  pyodideLoading = true;
  pyodideLoadPromise = (async () => {
    // Load Pyodide from CDN
    const script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/pyodide/v0.24.1/full/pyodide.js';
    document.head.appendChild(script);

    await new Promise((resolve) => {
      script.onload = resolve;
    });

    // Initialize Pyodide
    // eslint-disable-next-line no-undef
    pyodideInstance = await loadPyodide({
      indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.24.1/full/'
    });

    // Install NetworkX
    await pyodideInstance.loadPackage('networkx');

    pyodideLoading = false;
    return pyodideInstance;
  })();

  return pyodideLoadPromise;
}

/**
 * Uses NetworkX to check if graph is planar (without computing layout).
 * Returns { isPlanar: boolean, checked: boolean }
 */
async function checkPlanarity(nodes, edges) {
  // Very small graphs are always planar
  if (nodes.length < 3) {
    return { isPlanar: true, checked: true };
  }

  // Guard: skip heavy NetworkX work on very large graphs
  if (nodes.length > MAX_NODES_FOR_PLANARITY || edges.length > MAX_EDGES_FOR_PLANARITY) {
    console.warn(
      'Skipping NetworkX planarity check for large graph:',
      `nodes=${nodes.length}, edges=${edges.length}`
    );
    return { isPlanar: null, checked: false };
  }

  try {
    const pyodide = await loadPyodideAndNetworkX();

    const nodeIds = nodes.map(n => n.id);
    const edgeList = edges.map(e => [e.source, e.target]);

    const result = pyodide.runPython(`
import networkx as nx
import json

G = nx.Graph()
nodes = ${JSON.stringify(nodeIds)}
edges = ${JSON.stringify(edgeList)}

G.add_nodes_from(nodes)
G.add_edges_from(edges)

is_planar, _ = nx.check_planarity(G)
json.dumps({"isPlanar": is_planar})
`);

    const parsed = JSON.parse(result);
    return { isPlanar: parsed.isPlanar, checked: true };
  } catch (error) {
    console.error('Planarity check error:', error);
    return { isPlanar: null, checked: false };
  }
}

/**
 * Uses NetworkX to check if graph is planar and compute planar layout.
 * Returns { isPlanar, positions } where positions maps node id to {x, y}.
 */
async function computePlanarLayout(nodes, edges, width = 400, height = 400) {
  const pyodide = await loadPyodideAndNetworkX();

  // Convert graph to Python format
  const nodeIds = nodes.map(n => n.id);
  const edgeList = edges.map(e => [e.source, e.target]);

  // Run NetworkX planarity check and layout
  const result = pyodide.runPython(`
import networkx as nx
import json

# Create graph
G = nx.Graph()
nodes = ${JSON.stringify(nodeIds)}
edges = ${JSON.stringify(edgeList)}

G.add_nodes_from(nodes)
G.add_edges_from(edges)

# Check planarity and get embedding
is_planar, embedding = nx.check_planarity(G)

result = {"isPlanar": is_planar, "positions": None}

if is_planar:
    # Use the planar layout algorithm
    pos = nx.planar_layout(G, scale=1.0, center=(0, 0))
    # Convert numpy arrays to lists
    positions = {str(node): [float(coord[0]), float(coord[1])] for node, coord in pos.items()}
    result["positions"] = positions

json.dumps(result)
`);

  const parsed = JSON.parse(result);

  if (parsed.isPlanar && parsed.positions) {
    // Scale positions to fit our canvas
    const padding = 30;
    const positions = parsed.positions;
    
    // Find min/max to normalize
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const nodeId of nodeIds) {
      const pos = positions[String(nodeId)];
      if (pos) {
        minX = Math.min(minX, pos[0]);
        maxX = Math.max(maxX, pos[0]);
        minY = Math.min(minY, pos[1]);
        maxY = Math.max(maxY, pos[1]);
      }
    }

    const rangeX = maxX - minX || 1;
    const rangeY = maxY - minY || 1;

    // Map to canvas coordinates
    const scaledPositions = {};
    for (const nodeId of nodeIds) {
      const pos = positions[String(nodeId)];
      if (pos) {
        scaledPositions[nodeId] = {
          x: padding + ((pos[0] - minX) / rangeX) * (width - 2 * padding),
          y: padding + ((pos[1] - minY) / rangeY) * (height - 2 * padding)
        };
      }
    }

    return { isPlanar: true, positions: scaledPositions };
  }

  return { isPlanar: false, positions: null };
}

// ============================================================================
// GRID CONFIGURATION - Define the snap-to-grid system
// ============================================================================

const GRID_SIZE = 12;           // Size of each grid cell in pixels
const DEFAULT_GRID_COLS = 33;   // Default number of grid columns (33x33 = 1089 max nodes)
const DEFAULT_GRID_ROWS = 33;   // Default number of grid rows
const GRID_OFFSET = 8;          // Offset from edge to first grid point
// Grid points: 8, 20, 32, 44, ... 392 (33 points each axis, supports 1000+ nodes)

// Allow grid dimensions to adapt to imported file width/height
let GRID_COLS_DYNAMIC = DEFAULT_GRID_COLS;
let GRID_ROWS_DYNAMIC = DEFAULT_GRID_ROWS;

function getGridCols() {
  return GRID_COLS_DYNAMIC;
}

function getGridRows() {
  return GRID_ROWS_DYNAMIC;
}

/**
 * Updates dynamic grid dimensions based on an external canvas width/height.
 * This is used when importing a file that encodes its drawing area size.
 */
function updateGridFromSize(width, height) {
  if (!width || !height) return;

  const usableWidth = Math.max(0, width - 2 * GRID_OFFSET);
  const usableHeight = Math.max(0, height - 2 * GRID_OFFSET);

  const cols = Math.floor(usableWidth / GRID_SIZE) + 1;
  const rows = Math.floor(usableHeight / GRID_SIZE) + 1;

  // Keep within reasonable bounds to avoid performance issues
  GRID_COLS_DYNAMIC = Math.max(3, Math.min(cols, 100));
  GRID_ROWS_DYNAMIC = Math.max(3, Math.min(rows, 100));
}

/**
 * Snaps a coordinate to the nearest grid point.
 * @param {number} value - The coordinate to snap
 * @returns {number} - The snapped coordinate
 */
function snapToGrid(value) {
  const gridPoint = Math.round((value - GRID_OFFSET) / GRID_SIZE) * GRID_SIZE + GRID_OFFSET;
  const maxCoord = GRID_OFFSET + (getGridCols() - 1) * GRID_SIZE;
  return Math.max(GRID_OFFSET, Math.min(maxCoord, gridPoint));
}

/**
 * Gets all valid grid points as an array of {x, y} coordinates.
 * @returns {Array} - Array of grid point coordinates
 */
function getGridPoints() {
  const points = [];
  const rows = getGridRows();
  const cols = getGridCols();
  for (let row = 0; row < rows; row++) {
    for (let col = 0; col < cols; col++) {
      points.push({
        x: GRID_OFFSET + col * GRID_SIZE,
        y: GRID_OFFSET + row * GRID_SIZE
      });
    }
  }
  return points;
}

/**
 * Gets a random available grid point not already occupied by existing nodes.
 * @param {Array} existingNodes - Array of nodes already placed
 * @returns {Object|null} - {x, y} of available point, or null if none available
 */
function getRandomGridPoint(existingNodes) {
  const allPoints = getGridPoints();
  const availablePoints = allPoints.filter(point =>
    !existingNodes.some(node => node.x === point.x && node.y === point.y)
  );
  if (availablePoints.length === 0) return null;
  return availablePoints[Math.floor(Math.random() * availablePoints.length)];
}

/**
 * Snaps a node position to the nearest available grid point.
 * @param {number} x - X coordinate
 * @param {number} y - Y coordinate  
 * @param {Array} existingNodes - Other nodes to avoid
 * @param {number} currentNodeId - ID of current node (to exclude from collision check)
 * @returns {Object} - {x, y} snapped coordinates
 */
function snapToNearestAvailable(x, y, existingNodes, currentNodeId) {
  const snappedX = snapToGrid(x);
  const snappedY = snapToGrid(y);

  // Check if this grid point is occupied by another node
  const isOccupied = existingNodes.some(node =>
    node.id !== currentNodeId && node.x === snappedX && node.y === snappedY
  );

  if (!isOccupied) {
    return { x: snappedX, y: snappedY };
  }

  // Find nearest available grid point
  const allPoints = getGridPoints();
  let nearestPoint = { x: snappedX, y: snappedY };
  let nearestDist = Infinity;

  for (const point of allPoints) {
    const occupied = existingNodes.some(node =>
      node.id !== currentNodeId && node.x === point.x && node.y === point.y
    );
    if (!occupied) {
      const dist = Math.sqrt(Math.pow(point.x - x, 2) + Math.pow(point.y - y, 2));
      if (dist < nearestDist) {
        nearestDist = dist;
        nearestPoint = point;
      }
    }
  }

  return nearestPoint;
}

// ============================================================================
// SUBGRAPH DETECTION - Find K5 and K3,3 subgraphs
// ============================================================================

/**
 * Builds an adjacency list from edges for efficient neighbor lookup.
 * @param {Array} nodes - Array of node objects
 * @param {Array} edges - Array of edge objects with {source, target}
 * @returns {Map} - Map from node id to Set of neighbor ids
 */
function buildAdjacencyList(nodes, edges) {
  const adj = new Map();
  nodes.forEach(n => adj.set(n.id, new Set()));

  edges.forEach(e => {
    if (adj.has(e.source) && adj.has(e.target)) {
      adj.get(e.source).add(e.target);
      adj.get(e.target).add(e.source);
    }
  });

  return adj;
}

/**
 * Checks if a set of 5 nodes forms a K5 (complete graph on 5 vertices).
 * @param {Array} nodeIds - Array of 5 node IDs to check
 * @param {Map} adj - Adjacency list
 * @returns {boolean} - True if these 5 nodes form a K5
 */
function isK5(nodeIds, adj) {
  if (nodeIds.length !== 5) return false;

  // Check that every pair of nodes is connected
  for (let i = 0; i < 5; i++) {
    for (let j = i + 1; j < 5; j++) {
      if (!adj.get(nodeIds[i])?.has(nodeIds[j])) {
        return false;
      }
    }
  }
  return true;
}

/**
 * Checks if a set of 6 nodes forms a K3,3 (complete bipartite graph).
 * @param {Array} nodeIds - Array of 6 node IDs to check
 * @param {Map} adj - Adjacency list
 * @returns {Object|null} - {setA, setB} if K3,3 found, null otherwise
 */
function isK33(nodeIds, adj) {
  if (nodeIds.length !== 6) return null;

  // Try all possible ways to partition 6 nodes into two sets of 3
  // There are C(6,3) = 20 ways, but we only need to check 10 (half due to symmetry)
  const indices = [0, 1, 2, 3, 4, 5];

  for (let i = 0; i < 6; i++) {
    for (let j = i + 1; j < 6; j++) {
      for (let k = j + 1; k < 6; k++) {
        const setA = [nodeIds[i], nodeIds[j], nodeIds[k]];
        const setB = nodeIds.filter((_, idx) => idx !== i && idx !== j && idx !== k);

        // Check if every node in setA connects to every node in setB
        // and NO edges within setA or within setB
        let isValidK33 = true;

        // Check cross edges (must all exist)
        for (const a of setA) {
          for (const b of setB) {
            if (!adj.get(a)?.has(b)) {
              isValidK33 = false;
              break;
            }
          }
          if (!isValidK33) break;
        }

        if (isValidK33) {
          return { setA, setB };
        }
      }
    }
  }

  return null;
}

/**
 * Finds a K5 subgraph in the given graph.
 * Uses brute force for small graphs (checks all combinations of 5 nodes).
 * @param {Array} nodes - Array of node objects
 * @param {Array} edges - Array of edge objects
 * @returns {Object|null} - {nodes: [...], edges: [...]} of K5 subgraph, or null
 */
function findK5Subgraph(nodes, edges) {
  if (nodes.length < 5) return null;

  const adj = buildAdjacencyList(nodes, edges);
  const nodeIds = nodes.map(n => n.id);

  // For performance, limit search to first 20 nodes with highest degree
  const nodesByDegree = [...nodeIds].sort((a, b) =>
    (adj.get(b)?.size || 0) - (adj.get(a)?.size || 0)
  );
  const searchNodes = nodesByDegree.slice(0, Math.min(20, nodesByDegree.length));

  // Check all combinations of 5 nodes
  for (let i = 0; i < searchNodes.length; i++) {
    for (let j = i + 1; j < searchNodes.length; j++) {
      for (let k = j + 1; k < searchNodes.length; k++) {
        for (let l = k + 1; l < searchNodes.length; l++) {
          for (let m = l + 1; m < searchNodes.length; m++) {
            const candidates = [searchNodes[i], searchNodes[j], searchNodes[k], searchNodes[l], searchNodes[m]];
            if (isK5(candidates, adj)) {
              // Found K5! Return the nodes and edges
              const k5Edges = [];
              for (let a = 0; a < 5; a++) {
                for (let b = a + 1; b < 5; b++) {
                  k5Edges.push({ source: candidates[a], target: candidates[b] });
                }
              }
              return { nodes: candidates, edges: k5Edges };
            }
          }
        }
      }
    }
  }

  return null;
}

/**
 * Finds a K3,3 subgraph in the given graph.
 * @param {Array} nodes - Array of node objects
 * @param {Array} edges - Array of edge objects
 * @returns {Object|null} - {nodes: [...], edges: [...], setA, setB} of K3,3 subgraph, or null
 */
function findK33Subgraph(nodes, edges) {
  if (nodes.length < 6) return null;

  const adj = buildAdjacencyList(nodes, edges);
  const nodeIds = nodes.map(n => n.id);

  // For performance, limit search
  const searchNodes = nodeIds.slice(0, Math.min(15, nodeIds.length));

  // Check all combinations of 6 nodes
  for (let i = 0; i < searchNodes.length; i++) {
    for (let j = i + 1; j < searchNodes.length; j++) {
      for (let k = j + 1; k < searchNodes.length; k++) {
        for (let l = k + 1; l < searchNodes.length; l++) {
          for (let m = l + 1; m < searchNodes.length; m++) {
            for (let n = m + 1; n < searchNodes.length; n++) {
              const candidates = [searchNodes[i], searchNodes[j], searchNodes[k], searchNodes[l], searchNodes[m], searchNodes[n]];
              const result = isK33(candidates, adj);
              if (result) {
                // Found K3,3! Return the nodes and edges
                const k33Edges = [];
                for (const a of result.setA) {
                  for (const b of result.setB) {
                    k33Edges.push({ source: a, target: b });
                  }
                }
                return { nodes: candidates, edges: k33Edges, setA: result.setA, setB: result.setB };
              }
            }
          }
        }
      }
    }
  }

  return null;
}

/**
 * Finds a simple path between two nodes using BFS, avoiding certain nodes.
 * @param {number} start - Start node ID
 * @param {number} end - End node ID
 * @param {Map} adj - Adjacency list
 * @param {Set} avoidNodes - Set of node IDs to avoid (except start and end)
 * @param {number} maxLength - Maximum path length (default: 20)
 * @returns {Array|null} - Path as array of node IDs, or null if no path exists
 */
function findPath(start, end, adj, avoidNodes = new Set(), maxLength = 20) {
  if (start === end) return [start];
  
  const queue = [[start]];
  const visited = new Set([start]);
  
  while (queue.length > 0) {
    const path = queue.shift();
    const current = path[path.length - 1];
    
    // Limit path length to avoid very long paths
    if (path.length >= maxLength) continue;
    
    const neighbors = adj.get(current) || new Set();
    for (const neighbor of neighbors) {
      if (neighbor === end) {
        return [...path, neighbor];
      }
      
      if (!visited.has(neighbor) && !avoidNodes.has(neighbor)) {
        visited.add(neighbor);
        queue.push([...path, neighbor]);
      }
    }
  }
  
  return null;
}

/**
 * Finds a K5 subdivision in the given graph.
 * A K5 subdivision has 5 branch vertices with paths between each pair.
 * @param {Array} nodes - Array of node objects
 * @param {Array} edges - Array of edge objects
 * @returns {Object|null} - {nodes: [...], edges: [...], branchNodes: [...], paths: Map} or null
 */
function findK5Subdivision(nodes, edges) {
  if (nodes.length < 5) return null;
  
  const adj = buildAdjacencyList(nodes, edges);
  const nodeIds = nodes.map(n => n.id);
  
  // Try different sets of 5 nodes as potential branch vertices
  // Prioritize high-degree nodes as they're more likely to be branch vertices
  const nodesByDegree = [...nodeIds].sort((a, b) =>
    (adj.get(b)?.size || 0) - (adj.get(a)?.size || 0)
  );
  
  // Try up to 30 combinations of 5 nodes
  const maxCombinations = Math.min(30, Math.floor(nodesByDegree.length / 5));
  
  for (let attempt = 0; attempt < maxCombinations; attempt++) {
    // Select 5 candidate branch vertices
    const branchCandidates = [];
    const used = new Set();
    
    // Start with high-degree nodes
    for (const nodeId of nodesByDegree) {
      if (branchCandidates.length >= 5) break;
      if (!used.has(nodeId)) {
        branchCandidates.push(nodeId);
        used.add(nodeId);
      }
    }
    
    // If we don't have 5, fill with remaining nodes
    for (const nodeId of nodeIds) {
      if (branchCandidates.length >= 5) break;
      if (!used.has(nodeId)) {
        branchCandidates.push(nodeId);
        used.add(nodeId);
      }
    }
    
    if (branchCandidates.length < 5) break;
    
    // Try to find paths between all pairs of the 5 branch vertices
    const paths = new Map();
    const pathEdges = new Set();
    const allPathNodes = new Set(branchCandidates);
    let allPathsFound = true;
    
    // We need paths between all 10 pairs (5 choose 2)
    for (let i = 0; i < 5; i++) {
      for (let j = i + 1; j < 5; j++) {
        const start = branchCandidates[i];
        const end = branchCandidates[j];
        const key = `${Math.min(start, end)}-${Math.max(start, end)}`;
        
        // Avoid other branch vertices (except start and end)
        const avoidNodes = new Set(branchCandidates.filter(n => n !== start && n !== end));
        
        const path = findPath(start, end, adj, avoidNodes);
        
        if (!path || path.length < 2) {
          allPathsFound = false;
          break;
        }
        
        paths.set(key, path);
        
        // Collect all nodes and edges in the path
        for (let k = 0; k < path.length - 1; k++) {
          allPathNodes.add(path[k]);
          allPathNodes.add(path[k + 1]);
          pathEdges.add(`${Math.min(path[k], path[k + 1])}-${Math.max(path[k], path[k + 1])}`);
        }
      }
      if (!allPathsFound) break;
    }
    
    if (allPathsFound && paths.size === 10) {
      // Found a K5 subdivision!
      const subdivisionEdges = [];
      for (const edgeKey of pathEdges) {
        const [source, target] = edgeKey.split('-').map(Number);
        subdivisionEdges.push({ source, target });
      }
      
      return {
        nodes: Array.from(allPathNodes),
        edges: subdivisionEdges,
        branchNodes: branchCandidates,
        paths: paths,
        type: 'k5'
      };
    }
  }
  
  return null;
}

/**
 * Finds a K3,3 subdivision in the given graph.
 * A K3,3 subdivision has 6 branch vertices (3 in each partition) with paths between partitions.
 * @param {Array} nodes - Array of node objects
 * @param {Array} edges - Array of edge objects
 * @returns {Object|null} - {nodes: [...], edges: [...], branchNodes: {setA, setB}, paths: Map} or null
 */
function findK33Subdivision(nodes, edges) {
  if (nodes.length < 6) return null;
  
  const adj = buildAdjacencyList(nodes, edges);
  const nodeIds = nodes.map(n => n.id);
  
  // Try different sets of 6 nodes as potential branch vertices
  const nodesByDegree = [...nodeIds].sort((a, b) =>
    (adj.get(b)?.size || 0) - (adj.get(a)?.size || 0)
  );
  
  // Try up to 20 combinations
  const maxCombinations = Math.min(20, Math.floor(nodesByDegree.length / 6));
  
  for (let attempt = 0; attempt < maxCombinations; attempt++) {
    // Select 6 candidate branch vertices
    const branchCandidates = [];
    const used = new Set();
    
    for (const nodeId of nodesByDegree) {
      if (branchCandidates.length >= 6) break;
      if (!used.has(nodeId)) {
        branchCandidates.push(nodeId);
        used.add(nodeId);
      }
    }
    
    // Fill with remaining nodes if needed
    for (const nodeId of nodeIds) {
      if (branchCandidates.length >= 6) break;
      if (!used.has(nodeId)) {
        branchCandidates.push(nodeId);
        used.add(nodeId);
      }
    }
    
    if (branchCandidates.length < 6) break;
    
    // Try all possible ways to partition 6 nodes into two sets of 3
    for (let i = 0; i < 6; i++) {
      for (let j = i + 1; j < 6; j++) {
        for (let k = j + 1; k < 6; k++) {
          const setA = [branchCandidates[i], branchCandidates[j], branchCandidates[k]];
          const setB = branchCandidates.filter((_, idx) => idx !== i && idx !== j && idx !== k);
          
          // Try to find paths between all pairs (each node in setA to each node in setB)
          const paths = new Map();
          const pathEdges = new Set();
          const allPathNodes = new Set([...setA, ...setB]);
          let allPathsFound = true;
          
          // We need 9 paths (3 x 3)
          for (const a of setA) {
            for (const b of setB) {
              const key = `${Math.min(a, b)}-${Math.max(a, b)}`;
              
              // Avoid other branch vertices (except start and end)
              const avoidNodes = new Set(branchCandidates.filter(n => n !== a && n !== b));
              
              const path = findPath(a, b, adj, avoidNodes);
              
              if (!path || path.length < 2) {
                allPathsFound = false;
                break;
              }
              
              paths.set(key, path);
              
              // Collect all nodes and edges in the path
              for (let idx = 0; idx < path.length - 1; idx++) {
                allPathNodes.add(path[idx]);
                allPathNodes.add(path[idx + 1]);
                pathEdges.add(`${Math.min(path[idx], path[idx + 1])}-${Math.max(path[idx], path[idx + 1])}`);
              }
            }
            if (!allPathsFound) break;
          }
          
          if (allPathsFound && paths.size === 9) {
            // Found a K3,3 subdivision!
            const subdivisionEdges = [];
            for (const edgeKey of pathEdges) {
              const [source, target] = edgeKey.split('-').map(Number);
              subdivisionEdges.push({ source, target });
            }
            
            return {
              nodes: Array.from(allPathNodes),
              edges: subdivisionEdges,
              branchNodes: { setA, setB },
              paths: paths,
              type: 'k33'
            };
          }
        }
      }
    }
  }
  
  return null;
}

// ============================================================================
// UTILITY FUNCTIONS - Core algorithms for graph analysis
// ============================================================================

/**
 * Generates positions for nodes in a circular layout.
 * Used when uploaded graphs don't include position data.
 * 
 * @param {number} nodeCount - Number of nodes to position
 * @param {number} cx - Center X coordinate (default 200)
 * @param {number} cy - Center Y coordinate (default 200)
 * @param {number} radius - Radius of the circle (default 140)
 * @returns {Array} - Array of {x, y} positions
 */
function generateCircularLayout(nodeCount, cx = 200, cy = 200, radius = 140) {
  const positions = [];
  for (let i = 0; i < nodeCount; i++) {
    const angle = (i * 2 * Math.PI / nodeCount) - Math.PI / 2;
    positions.push({
      x: cx + radius * Math.cos(angle),
      y: cy + radius * Math.sin(angle)
    });
  }
  return positions;
}

/**
 * Parses a JSON graph file.
 * Supports multiple formats:
 * - { nodes: [{id, x, y}], edges: [{source, target}] }
 * - { nodes: [id, ...], edges: [[source, target], ...] }
 * - { vertices: [...], edges: [...] }
 * 
 * @param {string} content - JSON string content
 * @returns {Object} - { nodes, edges } or throws error
 */
function parseJSONGraph(content) {
  const data = JSON.parse(content);

  let nodes = [];
  let edges = [];

  // Handle nodes/vertices
  const rawNodes = data.nodes || data.vertices || data.Nodes || data.Vertices || [];

  if (rawNodes.length === 0) {
    throw new Error('No nodes found in file. Expected "nodes" or "vertices" array.');
  }

  // Check if nodes are objects with positions or just IDs
  if (typeof rawNodes[0] === 'object' && rawNodes[0] !== null) {
    // Nodes are objects - extract id and optional x, y
    nodes = rawNodes.map((n, index) => ({
      id: n.id !== undefined ? n.id : index,
      x: n.x !== undefined ? n.x : null,
      y: n.y !== undefined ? n.y : null
    }));
  } else {
    // Nodes are just IDs (numbers or strings)
    nodes = rawNodes.map((id, index) => ({
      id: typeof id === 'number' ? id : index,
      x: null,
      y: null
    }));
  }

  // Generate positions for nodes without coordinates
  const needsLayout = nodes.some(n => n.x === null || n.y === null);
  if (needsLayout) {
    const positions = generateCircularLayout(nodes.length);
    nodes = nodes.map((n, i) => ({
      ...n,
      x: n.x !== null ? n.x : positions[i].x,
      y: n.y !== null ? n.y : positions[i].y
    }));
  }

  // Scale positions to fit in our grid-aligned canvas
  const xs = nodes.map(n => n.x);
  const ys = nodes.map(n => n.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const rangeX = maxX - minX || 1;
  const rangeY = maxY - minY || 1;

  // Scale and snap to grid points. For performance on bigger uploads we avoid
  // the O(N^2) "nearest available" search and simply snap each node to the
  // nearest grid coordinate; mild overlaps are acceptable and can be fixed
  // interactively by the user.
  const cols = getGridCols();
  const rows = getGridRows();

  nodes = nodes.map(n => {
    const scaledX = GRID_OFFSET + ((n.x - minX) / rangeX) * (GRID_SIZE * (cols - 1));
    const scaledY = GRID_OFFSET + ((n.y - minY) / rangeY) * (GRID_SIZE * (rows - 1));
    return {
      ...n,
      x: snapToGrid(scaledX),
      y: snapToGrid(scaledY)
    };
  });

  // Create node ID mapping for edge resolution
  const nodeIdMap = new Map();
  nodes.forEach((n, index) => {
    nodeIdMap.set(n.id, index);
    nodeIdMap.set(String(n.id), index);
  });

  // Renumber nodes to sequential integers
  nodes = nodes.map((n, index) => ({ ...n, id: index }));

  // Handle edges
  const rawEdges = data.edges || data.links || data.Edges || data.Links || [];

  edges = rawEdges.map(e => {
    let source, target;

    if (Array.isArray(e)) {
      // Edge is [source, target]
      source = e[0];
      target = e[1];
    } else {
      // Edge is {source, target} or {from, to}
      source = e.source !== undefined ? e.source : e.from;
      target = e.target !== undefined ? e.target : e.to;
    }

    // Resolve to new node indices
    const sourceIdx = nodeIdMap.get(source) ?? nodeIdMap.get(String(source));
    const targetIdx = nodeIdMap.get(target) ?? nodeIdMap.get(String(target));

    if (sourceIdx === undefined || targetIdx === undefined) {
      return null; // Skip invalid edges
    }

    return { source: sourceIdx, target: targetIdx };
  }).filter(e => e !== null);

  return { nodes, edges };
}

/**
 * Parses a GraphML file.
 * 
 * @param {string} content - GraphML XML string
 * @returns {Object} - { nodes, edges } or throws error
 */
function parseGraphML(content) {
  const parser = new DOMParser();
  const doc = parser.parseFromString(content, 'text/xml');

  // Check for parse errors
  const parseError = doc.querySelector('parsererror');
  if (parseError) {
    throw new Error('Invalid GraphML: ' + parseError.textContent);
  }

  const graphEl = doc.querySelector('graph');
  if (!graphEl) {
    throw new Error('No <graph> element found in GraphML file.');
  }

  // Parse nodes
  const nodeEls = doc.querySelectorAll('node');
  if (nodeEls.length === 0) {
    throw new Error('No nodes found in GraphML file.');
  }

  const nodeIdMap = new Map();
  let nodes = [];

  nodeEls.forEach((nodeEl, index) => {
    const id = nodeEl.getAttribute('id') || String(index);
    nodeIdMap.set(id, index);

    // Try to get position from data elements (yFiles format)
    let x = null, y = null;
    const dataEls = nodeEl.querySelectorAll('data');
    dataEls.forEach(d => {
      const key = d.getAttribute('key');
      if (key === 'x' || key === 'd0') x = parseFloat(d.textContent);
      if (key === 'y' || key === 'd1') y = parseFloat(d.textContent);
    });

    // Also check for ShapeNode geometry (yFiles)
    const geometry = nodeEl.querySelector('Geometry, y\\:Geometry');
    if (geometry) {
      x = parseFloat(geometry.getAttribute('x')) || x;
      y = parseFloat(geometry.getAttribute('y')) || y;
    }

    nodes.push({ id: index, x, y });
  });

  // Generate positions for nodes without coordinates
  const needsLayout = nodes.some(n => n.x === null || n.y === null);
  if (needsLayout) {
    const positions = generateCircularLayout(nodes.length);
    nodes = nodes.map((n, i) => ({
      ...n,
      x: n.x !== null ? n.x : positions[i].x,
      y: n.y !== null ? n.y : positions[i].y
    }));
  }

  // Scale positions to fit canvas and snap to grid. As with JSON parsing,
  // we avoid the O(N^2) nearest-available search here for performance on
  // larger uploads and instead snap directly to the nearest grid point.
  const xs = nodes.map(n => n.x);
  const ys = nodes.map(n => n.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const rangeX = maxX - minX || 1;
  const rangeY = maxY - minY || 1;

  const cols = getGridCols();
  const rows = getGridRows();

  nodes = nodes.map(n => {
    const scaledX = GRID_OFFSET + ((n.x - minX) / rangeX) * (GRID_SIZE * (cols - 1));
    const scaledY = GRID_OFFSET + ((n.y - minY) / rangeY) * (GRID_SIZE * (rows - 1));
    return {
      ...n,
      x: snapToGrid(scaledX),
      y: snapToGrid(scaledY)
    };
  });

  // Parse edges
  const edgeEls = doc.querySelectorAll('edge');
  const edges = [];

  edgeEls.forEach(edgeEl => {
    const sourceId = edgeEl.getAttribute('source');
    const targetId = edgeEl.getAttribute('target');

    const sourceIdx = nodeIdMap.get(sourceId);
    const targetIdx = nodeIdMap.get(targetId);

    if (sourceIdx !== undefined && targetIdx !== undefined) {
      edges.push({ source: sourceIdx, target: targetIdx });
    }
  });

  return { nodes, edges };
}

/**
 * Parses a graph file based on its extension.
 * 
 * @param {File} file - The file object
 * @param {string} content - File content as string
 * @returns {Object} - { nodes, edges }
 */
function parseGraphFile(file, content) {
  const extension = file.name.split('.').pop().toLowerCase();

  if (extension === 'json') {
    return parseJSONGraph(content);
  } else if (extension === 'graphml' || extension === 'xml') {
    return parseGraphML(content);
  } else {
    // Try JSON first, then GraphML
    try {
      return parseJSONGraph(content);
    } catch {
      return parseGraphML(content);
    }
  }
}

/**
 * Checks if two line segments intersect using the cross product method.
 * Uses counter-clockwise (CCW) orientation to determine intersection.
 * 
 * @param {Object} p1 - Start point of first line segment (with id, x, y)
 * @param {Object} p2 - End point of first line segment
 * @param {Object} p3 - Start point of second line segment
 * @param {Object} p4 - End point of second line segment
 * @returns {boolean} - True if the segments intersect, false otherwise
 */
function doIntersect(p1, p2, p3, p4) {
  // If lines share a node (are neighbors), they don't count as crossing
  if (p1.id === p3.id || p1.id === p4.id || p2.id === p3.id || p2.id === p4.id) return false;

  // CCW helper: determines if three points are in counter-clockwise order
  const ccw = (a, b, c) => (c.y - a.y) * (b.x - a.x) > (b.y - a.y) * (c.x - a.x);

  // Two segments intersect if and only if the endpoints of each segment
  // are on opposite sides of the other segment
  return (ccw(p1, p3, p4) !== ccw(p2, p3, p4)) && (ccw(p1, p2, p3) !== ccw(p1, p2, p4));
}

/**
 * Calculates all edge crossings in a graph.
 * Compares every pair of edges to find intersections.
 * 
 * @param {Array} currentNodes - Array of node objects with {id, x, y}
 * @param {Array} currentEdges - Array of edge objects with {source, target}
 * @returns {Object} - Contains total crossings, max crossings per edge (k), and per-edge scores
 */
function calculateCrossings(currentNodes, currentEdges) {
  // Very small or empty graphs are trivial and cheap to analyze.
  let totalCrossings = 0;
  let scores = new Array(currentEdges.length).fill(0); // Track crossings per edge
  // Build a fast lookup map for node coordinates to avoid O(V) scans in the loop
  const nodeMap = new Map(currentNodes.map(n => [n.id, n]));

  // Compare edge pairs, but cap total comparisons to keep UI responsive on big graphs
  let comparisons = 0;

  const edgeCount = currentEdges.length;

  for (let i = 0; i < edgeCount; i++) {
    for (let j = i + 1; j < edgeCount; j++) {
      if (comparisons >= MAX_CROSSING_COMPARISONS) {
        // We’ve done enough work; return a partial but still useful picture
        const maxKPartial = scores.length > 0 ? Math.max(...scores) : 0;
        return { total: totalCrossings, k: maxKPartial, edgeScores: scores };
      }
      comparisons++;

      const edgeA = currentEdges[i];
      const edgeB = currentEdges[j];

      // Get the actual coordinates for each edge's endpoints
      const p1 = nodeMap.get(edgeA.source);
      const p2 = nodeMap.get(edgeA.target);
      const p3 = nodeMap.get(edgeB.source);
      const p4 = nodeMap.get(edgeB.target);

      // If edges intersect, increment counters
      if (p1 && p2 && p3 && p4 && doIntersect(p1, p2, p3, p4)) {
        totalCrossings++;
        scores[i]++; // Edge A has one more crossing
        scores[j]++; // Edge B has one more crossing
      }
    }
  }

  // k is the maximum number of crossings any single edge has
  const maxK = scores.length > 0 ? Math.max(...scores) : 0;
  return { total: totalCrossings, k: maxK, edgeScores: scores };
}

// ============================================================================
// GRAPH ANALYSIS HELPERS - For identifying specific graph types
// ============================================================================

/**
 * Checks if a graph is k-regular (all vertices have degree k).
 * @param {Array} nodes - Array of nodes
 * @param {Array} edges - Array of edges
 * @returns {number|boolean} - The degree k if regular, false otherwise
 */
function getRegularity(nodes, edges) {
  if (nodes.length === 0) return false;

  const adj = buildAdjacencyList(nodes, edges);
  const degrees = new Set();

  for (const neighbors of adj.values()) {
    degrees.add(neighbors.size);
  }

  // Also check nodes with no edges (might not be in adj map depending on implementation)
  if (adj.size < nodes.length) {
    degrees.add(0);
  }

  if (degrees.size === 1) {
    return degrees.values().next().value;
  }
  return false;
}

/**
 * Calculates the girth of a graph (length of the shortest cycle).
 * Uses BFS from each node to find the shortest cycle passing through it.
 * @param {Array} nodes - Array of nodes
 * @param {Array} edges - Array of edges
 * @returns {number} - The girth (Infinity if acyclic)
 */
function calculateGirth(nodes, edges) {
  const adj = buildAdjacencyList(nodes, edges);
  let minCycle = Infinity;

  for (const startNode of nodes) {
    const dist = new Map();
    const parent = new Map();
    const queue = [startNode.id];

    dist.set(startNode.id, 0);
    parent.set(startNode.id, null);

    while (queue.length > 0) {
      const u = queue.shift();

      // Optimization: if we've already gone deeper than current best/2, we can't find a shorter cycle
      if (dist.get(u) * 2 >= minCycle) break;

      const neighbors = adj.get(u) || new Set();
      for (const v of neighbors) {
        if (!dist.has(v)) {
          dist.set(v, dist.get(u) + 1);
          parent.set(v, u);
          queue.push(v);
        } else if (v !== parent.get(u)) {
          // Found a cycle
          const cycleLen = dist.get(u) + dist.get(v) + 1;
          minCycle = Math.min(minCycle, cycleLen);

          // Optimization: girth 3 is the minimum possible
          if (minCycle === 3) return 3;
        }
      }
    }
  }

  return minCycle;
}

/**
 * Identifies if the graph matches a known famous non-planar graph.
 * @param {Array} nodes - Array of nodes
 * @param {Array} edges - Array of edges
 * @returns {string|null} - The key of the famous graph (e.g., 'petersen') or null
 */
function identifyFamousGraph(nodes, edges) {
  const v = nodes.length;
  const e = edges.length;
  const k = getRegularity(nodes, edges);

  // Fast checks first (V, E, Regularity)

  // K5: 5 vertices, 10 edges, 4-regular
  if (v === 5 && e === 10 && k === 4) return 'k5';

  // K3,3: 6 vertices, 9 edges, 3-regular (and bipartite, but V/E/k is unique enough for small graphs)
  if (v === 6 && e === 9 && k === 3) {
    // Verify bipartite to be sure (no triangles)
    if (calculateGirth(nodes, edges) > 3) return 'k33';
  }

  // Petersen: 10 vertices, 15 edges, 3-regular, girth 5
  if (v === 10 && e === 15 && k === 3) {
    if (calculateGirth(nodes, edges) === 5) return 'petersen';
  }

  // Heawood: 14 vertices, 21 edges, 3-regular, girth 6
  if (v === 14 && e === 21 && k === 3) {
    if (calculateGirth(nodes, edges) === 6) return 'heawood';
  }

  // Kneser K(6,2): 15 vertices, 45 edges, 6-regular
  if (v === 15 && e === 45 && k === 6) return 'kneser62';

  // Mobius-Kantor: 16 vertices, 24 edges, 3-regular, girth 8
  if (v === 16 && e === 24 && k === 3) {
    if (calculateGirth(nodes, edges) === 8) return 'mobiusKantor';
  }

  // Pappus: 18 vertices, 27 edges, 3-regular, girth 6
  if (v === 18 && e === 27 && k === 3) {
    if (calculateGirth(nodes, edges) === 6) return 'pappus';
  }

  // Flower Snark J5: 20 vertices, 30 edges, 3-regular, girth 5 or 6?
  // Snarks have girth >= 5. J5 has girth 5.
  if (v === 20 && e === 30 && k === 3) {
    // Check girth to distinguish from other cubic graphs if needed, but for now V/E/k is strong hint
    const g = calculateGirth(nodes, edges);
    if (g >= 5) return 'flowerSnark';
  }

  return null;
}

/**
 * Finds nodes or edges whose removal would make the graph potentially planar.
 * "Potentially planar" means satisfying Euler's formula (E <= 3V-6) and having no K5/K3,3 subgraphs.
 * @param {Array} nodes - Current nodes
 * @param {Array} edges - Current edges
 * @returns {Object} - { nodes: Set<id>, edges: Set<id> } of valid fixes
 */
function findPlanarityFixes(nodes, edges) {
  const fixNodes = new Set();
  const fixEdges = new Set();

  // Helper to check planarity conditions
  const checkConditions = (testNodes, testEdges) => {
    const v = testNodes.length;
    const e = testEdges.length;
    if (v < 3) return true;
    if (e > 3 * v - 6) return false;
    if (findK5Subgraph(testNodes, testEdges)) return false;
    if (findK33Subgraph(testNodes, testEdges)) return false;
    return true;
  };

  // Try removing each edge
  // We need to identify edges uniquely. Using index in the array for now since that's stable during one check.
  // Actually, let's just store the edge objects or indices.
  // The UI will need to know which edge is which.
  // Let's store a signature "min-max" for edges.

  edges.forEach((edge, index) => {
    const remainingEdges = edges.filter((_, i) => i !== index);
    if (checkConditions(nodes, remainingEdges)) {
      fixEdges.add(`${Math.min(edge.source, edge.target)}-${Math.max(edge.source, edge.target)}`);
    }
  });

  // Try removing each node
  nodes.forEach(node => {
    const remainingNodes = nodes.filter(n => n.id !== node.id);
    const remainingEdges = edges.filter(e => e.source !== node.id && e.target !== node.id);
    if (checkConditions(remainingNodes, remainingEdges)) {
      fixNodes.add(node.id);
    }
  });

  return { nodes: fixNodes, edges: fixEdges };
}

/**
 * Generates all combinations of k elements from an array.
 * @param {Array} arr - The input array
 * @param {number} k - Number of elements per combination
 * @returns {Array} - Array of all combinations
 */
function combinations(arr, k) {
  if (k === 0) return [[]];
  if (k > arr.length) return [];
  
  const result = [];
  for (let i = 0; i <= arr.length - k; i++) {
    const head = arr[i];
    const tailCombos = combinations(arr.slice(i + 1), k - 1);
    for (const combo of tailCombos) {
      result.push([head, ...combo]);
    }
  }
  return result;
}

/**
 * Finds the minimum set of nodes whose removal makes the graph planar (verified with NetworkX).
 * @param {Array} nodes - Current nodes
 * @param {Array} edges - Current edges
 * @returns {Promise<{nodes: Array<number>, count: number}>} - Minimum set of node IDs to remove and count
 */
async function findMinimumNodeRemovalSet(nodes, edges) {
  // Skip expensive search on large graphs to keep the UI responsive
  if (nodes.length > MAX_NODES_FOR_MIN_REMOVAL || edges.length > MAX_EDGES_FOR_MIN_REMOVAL) {
    console.warn(
      'Skipping minimum node removal search for large graph:',
      `nodes=${nodes.length}, edges=${edges.length}`
    );
    return { nodes: [], count: Infinity };
  }

  const nodeIds = nodes.map(n => n.id);
  
  // Try removing 1 node, then 2, then 3, etc. until we find a solution
  for (let k = 1; k <= Math.min(nodeIds.length - 3, 5); k++) { // Limit to 5 nodes max for performance
    const combos = combinations(nodeIds, k);
    
    for (const nodeSet of combos) {
      const remainingNodes = nodes.filter(n => !nodeSet.includes(n.id));
      const remainingEdges = edges.filter(
        e => !nodeSet.includes(e.source) && !nodeSet.includes(e.target)
      );
      
      // Skip if graph becomes too small
      if (remainingNodes.length < 3) continue;
      
      // Check if the remaining graph is planar using NetworkX
      const result = await checkPlanarity(remainingNodes, remainingEdges);
      if (result.checked && result.isPlanar === true) {
        return { nodes: nodeSet, count: k };
      }
    }
  }
  
  return { nodes: [], count: Infinity };
}

/**
 * Finds the minimum set of edges whose removal makes the graph planar (verified with NetworkX).
 * @param {Array} nodes - Current nodes
 * @param {Array} edges - Current edges
 * @returns {Promise<{edges: Array<{source: number, target: number}>, count: number}>} - Minimum set of edges to remove and count
 */
async function findMinimumEdgeRemovalSet(nodes, edges) {
  // Skip expensive search on large graphs to keep the UI responsive
  if (nodes.length > MAX_NODES_FOR_MIN_REMOVAL || edges.length > MAX_EDGES_FOR_MIN_REMOVAL) {
    console.warn(
      'Skipping minimum edge removal search for large graph:',
      `nodes=${nodes.length}, edges=${edges.length}`
    );
    return { edges: [], count: Infinity };
  }

  // Try removing 1 edge, then 2, then 3, etc. until we find a solution
  for (let k = 1; k <= Math.min(edges.length, 10); k++) { // Limit to 10 edges max for performance
    const combos = combinations(edges, k);
    
    for (const edgeSet of combos) {
      // Create a set of edge keys for fast lookup
      const edgeSetKeys = new Set(edgeSet.map(e => 
        `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`
      ));
      
      const remainingEdges = edges.filter(e => {
        const key = `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`;
        return !edgeSetKeys.has(key);
      });
      
      // Skip if too few edges remain
      if (remainingEdges.length < 2) continue;
      
      // Check if the remaining graph is planar using NetworkX
      const result = await checkPlanarity(nodes, remainingEdges);
      if (result.checked && result.isPlanar === true) {
        return { edges: edgeSet, count: k };
      }
    }
  }
  
  return { edges: [], count: Infinity };
}

/**
 * Finds nodes whose removal actually makes the graph planar (verified with NetworkX).
 * @param {Array} nodes - Current nodes
 * @param {Array} edges - Current edges
 * @returns {Promise<Set<number>>} - Set of node IDs whose removal makes the graph planar
 */
async function findNodesThatMakePlanar(nodes, edges) {
  const planarNodes = new Set();
  
  // Check each node removal
  for (const node of nodes) {
    const remainingNodes = nodes.filter(n => n.id !== node.id);
    const remainingEdges = edges.filter(e => e.source !== node.id && e.target !== node.id);
    
    // Skip if graph becomes too small
    if (remainingNodes.length < 3) continue;
    
    // Check if the remaining graph is planar using NetworkX
    const result = await checkPlanarity(remainingNodes, remainingEdges);
    if (result.checked && result.isPlanar === true) {
      planarNodes.add(node.id);
    }
  }
  
  return planarNodes;
}

/**
 * Compares node removal vs edge removal and determines which is better.
 * @param {Object} nodeRemoval - {nodes: Array, count: number}
 * @param {Object} edgeRemoval - {edges: Array, count: number}
 * @returns {Object} - {better: 'nodes'|'edges'|'equal'|'none', reason: string}
 */
function compareRemovalStrategies(nodeRemoval, edgeRemoval) {
  // If neither works, return none
  if (nodeRemoval.count === Infinity && edgeRemoval.count === Infinity) {
    return { better: 'none', reason: 'Neither node nor edge removal can make this graph planar.' };
  }
  
  // If only one works, return that one
  if (nodeRemoval.count === Infinity) {
    return { 
      better: 'edges', 
      reason: `Removing ${edgeRemoval.count} edge${edgeRemoval.count > 1 ? 's' : ''} is the only viable option (removing nodes would require too many removals).` 
    };
  }
  
  if (edgeRemoval.count === Infinity) {
    return { 
      better: 'nodes', 
      reason: `Removing ${nodeRemoval.count} node${nodeRemoval.count > 1 ? 's' : ''} is the only viable option (removing edges would require too many removals).` 
    };
  }
  
  // Compare counts - fewer removals is better
  if (nodeRemoval.count < edgeRemoval.count) {
    return { 
      better: 'nodes', 
      reason: `Removing ${nodeRemoval.count} node${nodeRemoval.count > 1 ? 's' : ''} is better than removing ${edgeRemoval.count} edge${edgeRemoval.count > 1 ? 's' : ''} (fewer changes required).` 
    };
  }
  
  if (edgeRemoval.count < nodeRemoval.count) {
    return { 
      better: 'edges', 
      reason: `Removing ${edgeRemoval.count} edge${edgeRemoval.count > 1 ? 's' : ''} is better than removing ${nodeRemoval.count} node${nodeRemoval.count > 1 ? 's' : ''} (fewer changes required).` 
    };
  }
  
  // If equal, prefer edge removal (less destructive)
  return { 
    better: 'edges', 
    reason: `Both require ${nodeRemoval.count} removal${nodeRemoval.count > 1 ? 's' : ''}. Edge removal is preferred as it preserves all nodes.` 
  };
}

/**
 * Finds the best node to remove for planarity, selecting the node most involved in non-planarity.
 * @param {Array} nodes - Current nodes
 * @param {Array} edges - Current edges
 * @param {Set} fixNodes - Set of node IDs that can be removed to make graph planar
 * @param {Object|null} k5 - K5 subgraph info {nodes: [...], edges: [...]} or null
 * @param {Object|null} k33 - K3,3 subgraph info {nodes: [...], edges: [...], setA, setB} or null
 * @returns {number|null} - Node ID to remove, or null if no suitable node
 */
function findBestNodeToRemove(nodes, edges, fixNodes, k5, k33) {
  if (fixNodes.size === 0) return null;
  
  // Build adjacency list for degree calculation
  const adj = buildAdjacencyList(nodes, edges);
  
  // Count how many times each fixNode appears in K5/K3,3 subgraphs
  const involvementCount = new Map();
  const fixNodesArray = Array.from(fixNodes);
  
  // Initialize counts
  fixNodesArray.forEach(nodeId => {
    involvementCount.set(nodeId, 0);
  });
  
  // Count appearances in K5
  if (k5 && k5.nodes) {
    const k5NodeSet = new Set(k5.nodes);
    fixNodesArray.forEach(nodeId => {
      if (k5NodeSet.has(nodeId)) {
        involvementCount.set(nodeId, involvementCount.get(nodeId) + 1);
      }
    });
  }
  
  // Count appearances in K3,3
  if (k33 && k33.nodes) {
    const k33NodeSet = new Set(k33.nodes);
    fixNodesArray.forEach(nodeId => {
      if (k33NodeSet.has(nodeId)) {
        involvementCount.set(nodeId, involvementCount.get(nodeId) + 1);
      }
    });
  }
  
  // Find maximum involvement count
  let maxInvolvement = 0;
  fixNodesArray.forEach(nodeId => {
    const count = involvementCount.get(nodeId);
    if (count > maxInvolvement) {
      maxInvolvement = count;
    }
  });
  
  // Filter to nodes with maximum involvement
  const mostInvolved = fixNodesArray.filter(nodeId => 
    involvementCount.get(nodeId) === maxInvolvement
  );
  
  // If multiple nodes have same involvement, use highest degree as tiebreaker
  if (mostInvolved.length === 1) {
    return mostInvolved[0];
  }
  
  // Tiebreaker: highest degree
  let bestNode = mostInvolved[0];
  let maxDegree = adj.get(bestNode)?.size || 0;
  for (let i = 1; i < mostInvolved.length; i++) {
    const nodeId = mostInvolved[i];
    const degree = adj.get(nodeId)?.size || 0;
    if (degree > maxDegree) {
      maxDegree = degree;
      bestNode = nodeId;
    }
  }
  
  return bestNode;
}

/**
 * Generates an explanation for why specific nodes should be removed.
 * @param {Array} nodesToRemove - Node IDs to remove
 * @param {Array} allNodes - All nodes in the graph
 * @param {Array} allEdges - All edges in the graph
 * @param {Object|null} k5 - K5 subgraph info or null
 * @param {Object|null} k33 - K3,3 subgraph info or null
 * @returns {string} - Explanation text
 */
function generateRemovalExplanation(nodesToRemove, allNodes, allEdges, k5, k33) {
  const adj = buildAdjacencyList(allNodes, allEdges);
  const reasons = [];
  
  // Check involvement in K5
  if (k5 && k5.nodes) {
    const k5NodeSet = new Set(k5.nodes);
    const involvedInK5 = nodesToRemove.filter(id => k5NodeSet.has(id));
    if (involvedInK5.length > 0) {
      reasons.push(`${involvedInK5.length} node${involvedInK5.length > 1 ? 's' : ''} (${involvedInK5.join(', ')}) ${involvedInK5.length > 1 ? 'are' : 'is'} part of a K₅ subgraph`);
    }
  }
  
  // Check involvement in K3,3
  if (k33 && k33.nodes) {
    const k33NodeSet = new Set(k33.nodes);
    const involvedInK33 = nodesToRemove.filter(id => k33NodeSet.has(id));
    if (involvedInK33.length > 0) {
      reasons.push(`${involvedInK33.length} node${involvedInK33.length > 1 ? 's' : ''} (${involvedInK33.join(', ')}) ${involvedInK33.length > 1 ? 'are' : 'is'} part of a K₃,₃ subgraph`);
    }
  }
  
  // Check degrees (high degree nodes are often problematic)
  const degrees = nodesToRemove.map(id => ({ id, degree: adj.get(id)?.size || 0 }));
  degrees.sort((a, b) => b.degree - a.degree);
  if (degrees.length > 0 && degrees[0].degree > 3) {
    const highDegreeNodes = degrees.filter(d => d.degree > 3).map(d => d.id);
    if (highDegreeNodes.length > 0) {
      reasons.push(`${highDegreeNodes.length} node${highDegreeNodes.length > 1 ? 's' : ''} (${highDegreeNodes.join(', ')}) ${highDegreeNodes.length > 1 ? 'have' : 'has'} high degree (${degrees[0].degree} connections)`);
    }
  }
  
  if (reasons.length === 0) {
    return `Removing these ${nodesToRemove.length} node${nodesToRemove.length > 1 ? 's' : ''} (${nodesToRemove.join(', ')}) will break the non-planar structure.`;
  }
  
  return `These nodes are being removed because: ${reasons.join('; ')}.`;
}

/**
 * Generates an explanation for why specific edges should be removed.
 * @param {Array} edgesToRemove - Edges to remove [{source, target}, ...]
 * @param {Array} allNodes - All nodes in the graph
 * @param {Array} allEdges - All edges in the graph
 * @param {Object|null} k5 - K5 subgraph info or null
 * @param {Object|null} k33 - K3,3 subgraph info or null
 * @returns {string} - Explanation text
 */
function generateEdgeRemovalExplanation(edgesToRemove, allNodes, allEdges, k5, k33) {
  const reasons = [];
  const edgeStrings = edgesToRemove.map(e => `${e.source}-${e.target}`);
  
  // Check if edges are part of K5
  if (k5 && k5.edges) {
    const k5EdgeSet = new Set(k5.edges.map(e => `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`));
    const involvedInK5 = edgesToRemove.filter(e => {
      const key = `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`;
      return k5EdgeSet.has(key);
    });
    if (involvedInK5.length > 0) {
      reasons.push(`${involvedInK5.length} edge${involvedInK5.length > 1 ? 's' : ''} ${involvedInK5.length > 1 ? 'are' : 'is'} part of a K₅ subgraph`);
    }
  }
  
  // Check if edges are part of K3,3
  if (k33 && k33.edges) {
    const k33EdgeSet = new Set(k33.edges.map(e => `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`));
    const involvedInK33 = edgesToRemove.filter(e => {
      const key = `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`;
      return k33EdgeSet.has(key);
    });
    if (involvedInK33.length > 0) {
      reasons.push(`${involvedInK33.length} edge${involvedInK33.length > 1 ? 's' : ''} ${involvedInK33.length > 1 ? 'are' : 'is'} part of a K₃,₃ subgraph`);
    }
  }
  
  if (reasons.length === 0) {
    return `Removing these ${edgesToRemove.length} edge${edgesToRemove.length > 1 ? 's' : ''} (${edgeStrings.join(', ')}) will break the non-planar structure.`;
  }
  
  return `These edges are being removed because: ${reasons.join('; ')}.`;
}

// ============================================================================
// EXAMPLE GRAPHS - Famous non-planar graphs for demonstration
// ============================================================================

/**
 * Predefined example graphs that are known to be non-planar.
 * Each graph includes metadata and a generator function to create its structure.
 */
const EXAMPLE_GRAPHS = {
  // K₅: Complete graph on 5 vertices (every vertex connected to every other)
  k5: {
    name: 'K₅ (Complete Graph)',
    description: 'The complete graph on 5 vertices',
    explanation: 'K₅ is non-planar. By Kuratowski\'s theorem, any graph containing K₅ as a subgraph cannot be drawn on a plane without edge crossings. This is one of the two fundamental non-planar graphs.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      // Place 5 nodes on grid points in a pentagon-like formation
      // Using snapToGrid to get valid coordinates
      const positions = [
        { x: snapToGrid(200), y: snapToGrid(44) },   // Top center
        { x: snapToGrid(320), y: snapToGrid(140) },  // Top right
        { x: snapToGrid(284), y: snapToGrid(308) },  // Bottom right
        { x: snapToGrid(116), y: snapToGrid(308) },  // Bottom left
        { x: snapToGrid(80), y: snapToGrid(140) },   // Top left
      ];

      for (let i = 0; i < 5; i++) {
        nodes.push({ id: i, x: positions[i].x, y: positions[i].y });
      }

      // Connect every pair of nodes (complete graph)
      for (let i = 0; i < 5; i++) {
        for (let j = i + 1; j < 5; j++) {
          edges.push({ source: i, target: j });
        }
      }
      return { nodes, edges };
    }
  },

  // K₃,₃: Complete bipartite graph (two sets of 3 vertices, all cross-connections)
  k33: {
    name: 'K₃,₃ (Bipartite)',
    description: 'Complete bipartite graph',
    explanation: 'K₃,₃ is non-planar. It\'s the complete bipartite graph with two sets of 3 vertices each, where every vertex in one set connects to all vertices in the other. Along with K₅, it forms the basis of Kuratowski\'s theorem.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      // Create two columns of 3 nodes each on grid points
      const leftX = snapToGrid(116);
      const rightX = snapToGrid(284);
      const yPositions = [snapToGrid(80), snapToGrid(200), snapToGrid(320)];

      for (let i = 0; i < 3; i++) {
        nodes.push({ id: i, x: leftX, y: yPositions[i] });         // Left column
        nodes.push({ id: i + 3, x: rightX, y: yPositions[i] });    // Right column
      }

      // Connect every node in left column to every node in right column
      for (let i = 0; i < 3; i++) {
        for (let j = 3; j < 6; j++) {
          edges.push({ source: i, target: j });
        }
      }
      return { nodes, edges };
    }
  },

  // Petersen Graph: Famous 10-vertex graph used as counterexample in graph theory
  petersen: {
    name: 'Petersen Graph',
    description: 'Famous non-planar graph',
    explanation: 'The Petersen graph is non-planar. It\'s a famous graph in mathematics with 10 vertices and 15 edges. It contains K₅ as a minor, making it non-planar. It\'s often used as a counterexample in graph theory.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      // Outer pentagon on grid points
      const outerPositions = [
        { x: snapToGrid(200), y: snapToGrid(32) },   // Top
        { x: snapToGrid(332), y: snapToGrid(128) },  // Top right
        { x: snapToGrid(284), y: snapToGrid(320) },  // Bottom right
        { x: snapToGrid(116), y: snapToGrid(320) },  // Bottom left
        { x: snapToGrid(68), y: snapToGrid(128) },   // Top left
      ];

      // Inner pentagon on grid points
      const innerPositions = [
        { x: snapToGrid(200), y: snapToGrid(116) },  // Top
        { x: snapToGrid(272), y: snapToGrid(164) },  // Top right
        { x: snapToGrid(248), y: snapToGrid(260) },  // Bottom right
        { x: snapToGrid(152), y: snapToGrid(260) },  // Bottom left
        { x: snapToGrid(128), y: snapToGrid(164) },  // Top left
      ];

      // Add outer nodes (0-4)
      for (let i = 0; i < 5; i++) {
        nodes.push({ id: i, x: outerPositions[i].x, y: outerPositions[i].y });
      }

      // Add inner nodes (5-9)
      for (let i = 0; i < 5; i++) {
        nodes.push({ id: i + 5, x: innerPositions[i].x, y: innerPositions[i].y });
      }

      // Create the Petersen graph's unique edge structure
      for (let i = 0; i < 5; i++) {
        edges.push({ source: i, target: (i + 1) % 5 });           // Outer pentagon edges
        edges.push({ source: i, target: i + 5 });                  // Spokes connecting to inner
        edges.push({ source: i + 5, target: ((i + 2) % 5) + 5 }); // Inner pentagram (skip 1)
      }
      return { nodes, edges };
    }
  },

  // Flower Snark J5: A famous snark graph (cyclically 4-edge-connected cubic graph with chromatic index 4)
  flowerSnark: {
    name: 'Flower Snark J₅',
    description: '20-vertex snark graph',
    explanation: 'The Flower Snark J₅ is a non-planar snark with 20 vertices and 30 edges. Snarks are bridgeless cubic graphs that require 4 colors to edge-color. This graph contains K₃,₃ as a minor, making it non-planar.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      // Create 4 groups of 5 nodes each in a flower-like pattern
      const cx = 200, cy = 200;
      const outerRadius = 160;
      const innerRadius = 80;

      // Outer ring (0-4)
      for (let i = 0; i < 5; i++) {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2;
        nodes.push({
          id: i,
          x: snapToGrid(cx + outerRadius * Math.cos(angle)),
          y: snapToGrid(cy + outerRadius * Math.sin(angle))
        });
      }

      // Middle ring A (5-9)
      for (let i = 0; i < 5; i++) {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2 + Math.PI / 10;
        nodes.push({
          id: i + 5,
          x: snapToGrid(cx + (outerRadius * 0.7) * Math.cos(angle)),
          y: snapToGrid(cy + (outerRadius * 0.7) * Math.sin(angle))
        });
      }

      // Middle ring B (10-14)
      for (let i = 0; i < 5; i++) {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2 - Math.PI / 10;
        nodes.push({
          id: i + 10,
          x: snapToGrid(cx + (innerRadius * 1.2) * Math.cos(angle)),
          y: snapToGrid(cy + (innerRadius * 1.2) * Math.sin(angle))
        });
      }

      // Inner ring (15-19)
      for (let i = 0; i < 5; i++) {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2;
        nodes.push({
          id: i + 15,
          x: snapToGrid(cx + innerRadius * 0.5 * Math.cos(angle)),
          y: snapToGrid(cy + innerRadius * 0.5 * Math.sin(angle))
        });
      }

      // Flower snark edges
      for (let i = 0; i < 5; i++) {
        // Outer ring connections
        edges.push({ source: i, target: i + 5 });
        edges.push({ source: i, target: i + 10 });
        // Middle connections
        edges.push({ source: i + 5, target: (i + 1) % 5 + 5 });
        edges.push({ source: i + 5, target: i + 15 });
        edges.push({ source: i + 10, target: (i + 4) % 5 + 10 });
        edges.push({ source: i + 10, target: i + 15 });
      }

      return { nodes, edges };
    }
  },

  // Heawood Graph: The (3,6)-cage, smallest graph with girth 6 that is 3-regular
  heawood: {
    name: 'Heawood Graph',
    description: '(3,6)-cage with 14 vertices',
    explanation: 'The Heawood graph is the (3,6)-cage: the smallest 3-regular graph with girth 6. It has 14 vertices and 21 edges. It\'s the incidence graph of the Fano plane and is non-planar, containing K₃,₃ as a subdivision.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      const cx = 200, cy = 200;
      const radius = 150;

      // All 14 vertices arranged in a circle
      for (let i = 0; i < 14; i++) {
        const angle = (i * 2 * Math.PI / 14) - Math.PI / 2;
        nodes.push({
          id: i,
          x: snapToGrid(cx + radius * Math.cos(angle)),
          y: snapToGrid(cy + radius * Math.sin(angle))
        });
      }

      // Heawood graph edges:
      // 1. The outer 14-cycle
      for (let i = 0; i < 14; i++) {
        edges.push({ source: i, target: (i + 1) % 14 });
      }

      // 2. The 7 chord edges that create crossings
      // Connect alternating vertices with skip-5 pattern
      for (let i = 0; i < 7; i++) {
        edges.push({ source: i * 2, target: (i * 2 + 5) % 14 });
      }
      // Total: 14 + 7 = 21 edges, each vertex has degree 3

      return { nodes, edges };
    }
  },

  // Möbius-Kantor Graph: The (3,8)-cage
  mobiusKantor: {
    name: 'Möbius-Kantor Graph',
    description: '(3,8)-cage with 16 vertices',
    explanation: 'The Möbius-Kantor graph is the (3,8)-cage: the smallest 3-regular graph with girth 8. It has 16 vertices and 24 edges. It\'s a generalized Petersen graph GP(8,3) and is non-planar.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      const cx = 200, cy = 200;
      const outerRadius = 150;
      const innerRadius = 70;

      // Outer octagon (8 vertices)
      for (let i = 0; i < 8; i++) {
        const angle = (i * 2 * Math.PI / 8) - Math.PI / 2;
        nodes.push({
          id: i,
          x: snapToGrid(cx + outerRadius * Math.cos(angle)),
          y: snapToGrid(cy + outerRadius * Math.sin(angle))
        });
      }

      // Inner octagon (8 vertices)
      for (let i = 0; i < 8; i++) {
        const angle = (i * 2 * Math.PI / 8) - Math.PI / 2;
        nodes.push({
          id: i + 8,
          x: snapToGrid(cx + innerRadius * Math.cos(angle)),
          y: snapToGrid(cy + innerRadius * Math.sin(angle))
        });
      }

      // Möbius-Kantor edges (generalized Petersen GP(8,3))
      for (let i = 0; i < 8; i++) {
        edges.push({ source: i, target: (i + 1) % 8 });              // Outer cycle
        edges.push({ source: i, target: i + 8 });                     // Spokes
        edges.push({ source: i + 8, target: ((i + 3) % 8) + 8 });    // Inner star (step 3)
      }

      return { nodes, edges };
    }
  },

  // Kneser Graph K(6,2): 15 vertices representing 2-subsets of {1,2,3,4,5,6}
  kneser62: {
    name: 'Kneser K(6,2)',
    description: '15-vertex Kneser graph',
    explanation: 'The Kneser graph K(6,2) has 15 vertices representing all 2-element subsets of a 6-element set. Two vertices are connected if their subsets are disjoint. It\'s non-planar and highly symmetric.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      // Generate all 2-subsets of {0,1,2,3,4,5}
      const subsets = [];
      for (let i = 0; i < 6; i++) {
        for (let j = i + 1; j < 6; j++) {
          subsets.push([i, j]);
        }
      }

      const cx = 200, cy = 200;

      // Arrange in three rings
      const ring1 = [0, 1, 2, 3, 4];       // 5 vertices outer
      const ring2 = [5, 6, 7, 8, 9];       // 5 vertices middle
      const ring3 = [10, 11, 12, 13, 14];  // 5 vertices inner

      const radii = [150, 100, 50];

      ring1.forEach((idx, i) => {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2;
        nodes.push({
          id: idx,
          x: snapToGrid(cx + radii[0] * Math.cos(angle)),
          y: snapToGrid(cy + radii[0] * Math.sin(angle))
        });
      });

      ring2.forEach((idx, i) => {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2 + Math.PI / 5;
        nodes.push({
          id: idx,
          x: snapToGrid(cx + radii[1] * Math.cos(angle)),
          y: snapToGrid(cy + radii[1] * Math.sin(angle))
        });
      });

      ring3.forEach((idx, i) => {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2;
        nodes.push({
          id: idx,
          x: snapToGrid(cx + radii[2] * Math.cos(angle)),
          y: snapToGrid(cy + radii[2] * Math.sin(angle))
        });
      });

      // Connect disjoint subsets
      for (let i = 0; i < subsets.length; i++) {
        for (let j = i + 1; j < subsets.length; j++) {
          const [a1, a2] = subsets[i];
          const [b1, b2] = subsets[j];
          // Check if subsets are disjoint
          if (a1 !== b1 && a1 !== b2 && a2 !== b1 && a2 !== b2) {
            edges.push({ source: i, target: j });
          }
        }
      }

      return { nodes, edges };
    }
  },

  // Pappus Graph: A famous (3,6)-bipartite cage
  pappus: {
    name: 'Pappus Graph',
    description: '18-vertex symmetric graph',
    explanation: 'The Pappus graph is a bipartite 3-regular graph with 18 vertices and 27 edges. It\'s named after Pappus of Alexandria and is related to the Pappus configuration in projective geometry. It\'s non-planar.',
    isPlanar: false,
    generate: () => {
      const nodes = [];
      const edges = [];

      const cx = 200, cy = 200;
      const outerRadius = 160;
      const middleRadius = 100;
      const innerRadius = 45;

      // Three hexagons
      for (let ring = 0; ring < 3; ring++) {
        const radius = ring === 0 ? outerRadius : ring === 1 ? middleRadius : innerRadius;
        const offset = ring * Math.PI / 6;
        for (let i = 0; i < 6; i++) {
          const angle = (i * 2 * Math.PI / 6) - Math.PI / 2 + offset;
          nodes.push({
            id: ring * 6 + i,
            x: snapToGrid(cx + radius * Math.cos(angle)),
            y: snapToGrid(cy + radius * Math.sin(angle))
          });
        }
      }

      // Pappus graph edges
      // Outer hexagon cycle
      for (let i = 0; i < 6; i++) {
        edges.push({ source: i, target: (i + 1) % 6 });
      }
      // Inner hexagon cycle  
      for (let i = 0; i < 6; i++) {
        edges.push({ source: 12 + i, target: 12 + (i + 1) % 6 });
      }
      // Connections between rings
      for (let i = 0; i < 6; i++) {
        edges.push({ source: i, target: 6 + i });
        edges.push({ source: 6 + i, target: 12 + i });
        edges.push({ source: 6 + i, target: 12 + (i + 3) % 6 });
      }

      return { nodes, edges };
    }
  }
};

// ============================================================================
// MAIN APP COMPONENT
// ============================================================================

function App() {
  // --------------------------------------------------------------------------
  // STATE MANAGEMENT
  // --------------------------------------------------------------------------

  // Graph data state
  const [nodes, setNodes] = useState([]);                    // Array of {id, x, y} objects
  const [edges, setEdges] = useState([]);                    // Array of {source, target} objects
  const [nodeCountInput, setNodeCountInput] = useState('20'); // User input for node count
  const [metrics, setMetrics] = useState({ total: 0, k: 0, edgeScores: [] }); // Crossing analysis
  const [networkxPlanarity, setNetworkxPlanarity] = useState({ isPlanar: null, checked: false, checking: false }); // NetworkX result

  // UI state
  const [isDragOver, setIsDragOver] = useState(false);       // File drag-drop hover state
  const [showGlobalLoading, setShowGlobalLoading] = useState(false); // Full-app loading overlay
  const [activeAnalysis, setActiveAnalysis] = useState(null); // Currently displayed analysis panel
  const [currentChallenge, setCurrentChallenge] = useState(null); // Active challenge graph ID
  const [showSolution, setShowSolution] = useState(false);   // Solution visibility toggle
  const [uploadStatus, setUploadStatus] = useState(null);    // { type: 'success'|'error', message: string }
  const [uploadedFileName, setUploadedFileName] = useState(null); // Name of uploaded file
  const [isUploading, setIsUploading] = useState(false);     // Prevents concurrent uploads and shows loading state
  const [highlightedSubgraph, setHighlightedSubgraph] = useState(null); // 'k5' | 'k33' | null
  const [showFixes, setShowFixes] = useState(false); // Toggle for fix suggestions
  const [originalPositions, setOriginalPositions] = useState(null); // Store original node positions before animation
  const [isAnimating, setIsAnimating] = useState(false); // Animation in progress
  const [currentlyMovingNode, setCurrentlyMovingNode] = useState(null); // Node ID currently being animated
  const [show3DView, setShow3DView] = useState(false); // 3D sphere view toggle
  const [sphereRotation, setSphereRotation] = useState({ x: 0, y: 0 }); // 3D rotation angles
  const [isDragging3D, setIsDragging3D] = useState(false); // 3D view drag state
  const [autoRotate3D, setAutoRotate3D] = useState(true); // Auto-rotate toggle
  const [lastMousePos, setLastMousePos] = useState({ x: 0, y: 0 }); // For drag delta calculation
  const [loadedExampleId, setLoadedExampleId] = useState(null); // Track which example graph was loaded
  const [isPanning2D, setIsPanning2D] = useState(false); // 2D graph pan state
  const [panStart, setPanStart] = useState({ x: 0, y: 0 }); // Pan start position
  const [planarityCheckTime, setPlanarityCheckTime] = useState(null); // Time taken for planarity check in ms
  const [showThoughtProcess, setShowThoughtProcess] = useState(false); // Toggle for thought process display

  // Minimum removal sets for planarity
  const [minimumRemovalSet, setMinimumRemovalSet] = useState(null); // {nodes: Array, count: number, explanation: string} or null
  const [minimumEdgeRemovalSet, setMinimumEdgeRemovalSet] = useState(null); // {edges: Array, count: number, explanation: string} or null
  const [removalComparison, setRemovalComparison] = useState(null); // {better: 'nodes'|'edges'|'equal'|'none', reason: string} or null
  const [checkingMinimumSet, setCheckingMinimumSet] = useState(false); // Whether we're currently computing minimum sets

  // Zoom and pan state
  // Start slightly zoomed out at 80% (showing more of the grid)
  const [zoom, setZoom] = useState(0.8);                     // Zoom level (1 = 100%)
  const [panOffset, setPanOffset] = useState({ x: 0, y: 0 }); // Pan offset for panning the view

  // Grid and UI preferences
  const [showGrid, setShowGrid] = useState(true);            // Show/hide grid
  const [gridStyle, setGridStyle] = useState('dots');        // 'dots', 'lines', 'both'
  const [hoveredEdge, setHoveredEdge] = useState(null);     // Currently hovered edge index
  const [selectedNode, setSelectedNode] = useState(null);   // Selected node ID
  const [showExportMenu, setShowExportMenu] = useState(false); // Export menu visibility
  const [isDraggingNode, setIsDraggingNode] = useState(false); // Whether a node is currently being dragged

  // --------------------------------------------------------------------------
  // REFS - For DOM access and mutable values that don't trigger re-renders
  // --------------------------------------------------------------------------

  const svgRef = useRef(null);           // Reference to the interactive SVG element
  const fileInputRef = useRef(null);     // Reference to hidden file input
  const draggingNodeRef = useRef(null);  // Currently dragged node ID (using ref to avoid stale closures)
  const dragOffsetRef = useRef({ x: 0, y: 0 }); // Offset between mouse and node center when dragging starts
  const dragStartPosRef = useRef({ x: 0, y: 0 }); // Initial mouse position when mousedown occurs
  const isDraggingRef = useRef(false);  // Whether actual dragging has started (after threshold)
  const edgesRef = useRef(edges);        // Ref to current edges for use in callbacks

  // Keep edgesRef in sync with edges state
  edgesRef.current = edges;

  // --------------------------------------------------------------------------
  // EXAMPLE GRAPHS LIST - For the cards in the UI
  // --------------------------------------------------------------------------

  const [importedGraphs] = useState([
    { id: 'k5', name: 'K₅ (Complete Graph)', nodes: 5, edges: 10, description: 'The complete graph on 5 vertices' },
    { id: 'k33', name: 'K₃,₃ (Bipartite)', nodes: 6, edges: 9, description: 'Complete bipartite graph' },
    { id: 'petersen', name: 'Petersen Graph', nodes: 10, edges: 15, description: 'Famous non-planar graph' },
    { id: 'flowerSnark', name: 'Flower Snark J₅', nodes: 20, edges: 30, description: '20-vertex snark graph' },
    { id: 'heawood', name: 'Heawood Graph', nodes: 14, edges: 21, description: '(3,6)-cage with 14 vertices' },
    { id: 'mobiusKantor', name: 'Möbius-Kantor', nodes: 16, edges: 24, description: '(3,8)-cage with 16 vertices' },
    { id: 'kneser62', name: 'Kneser K(6,2)', nodes: 15, edges: 45, description: '15-vertex Kneser graph' },
    { id: 'pappus', name: 'Pappus Graph', nodes: 18, edges: 27, description: '18-vertex symmetric graph' },
  ]);

  // --------------------------------------------------------------------------
  // GRAPH ANALYSIS FUNCTIONS
  // --------------------------------------------------------------------------

  /**
   * Analyzes the current graph and updates metrics state.
   * Wrapped in useCallback to maintain referential stability.
   */
  const analyzeGraph = useCallback((currentNodes, currentEdges) => {
    const result = calculateCrossings(currentNodes, currentEdges);
    setMetrics(result);
  }, []);

  /**
   * Generates a random graph with the specified number of nodes.
   * Nodes are placed on grid points. Each node gets 1-2 random edges.
   */
  const generateRandomGraph = useCallback((nodeCount = 10) => {
    const maxNodes = getGridCols() * getGridRows(); // Maximum nodes that can fit on grid
    const actualNodeCount = Math.min(nodeCount, maxNodes);

    const newNodes = [];
    const newEdges = [];

    // Create nodes at random grid positions
    for (let i = 0; i < actualNodeCount; i++) {
      const point = getRandomGridPoint(newNodes);
      if (point) {
        newNodes.push({
          id: i,
          x: point.x,
          y: point.y,
        });
      }
    }

    // Create random edges (each node connects to 1-2 other random nodes)
    for (let i = 0; i < newNodes.length; i++) {
      const target1 = Math.floor(Math.random() * newNodes.length);
      const target2 = Math.floor(Math.random() * newNodes.length);
      if (i !== target1) newEdges.push({ source: i, target: target1 });
      if (i !== target2 && target1 !== target2) newEdges.push({ source: i, target: target2 });
    }

    setNodes(newNodes);
    setEdges(newEdges);
    setLoadedExampleId(null);  // Clear any loaded example
    setActiveAnalysis(null);  // Close example graph analysis panel
    setCurrentChallenge(null);
    setShowFixes(false); // Reset fixes visibility
    analyzeGraph(newNodes, newEdges);
  }, [analyzeGraph]);

  /**
   * Generates a fully random graph with random node count (5-15 nodes).
   */
  const generateFullyRandom = () => {
    const randomNodeCount = Math.floor(Math.random() * 11) + 5;
    setNodeCountInput(String(randomNodeCount));
    generateRandomGraph(randomNodeCount);
  };

  /**
   * Generates a graph with the user-specified node count.
   * Validates input and enforces min/max limits.
   */
  const generateWithNodeCount = () => {
    const count = parseInt(nodeCountInput, 10);
    if (isNaN(count) || count < 2) return;  // Minimum 2 nodes
    if (count > 1000) return;                // Maximum 1000 nodes
    generateRandomGraph(count);
  };

  // Initial state is empty - user must upload or generate a graph

  // Check planarity using NetworkX whenever graph changes
  useEffect(() => {
    if (nodes.length < 2 || edges.length === 0) {
      setNetworkxPlanarity({ isPlanar: true, checked: true, checking: false });
      setPlanarityCheckTime(null);
      return;
    }

    // Skip if we're currently dragging

    // Create a stable key for this graph configuration
    const graphKey = `${nodes.length}-${edges.length}-${edges.map(e => `${e.source}-${e.target}`).join(',')}`;
    
    setNetworkxPlanarity(prev => ({ ...prev, checking: true }));

    // Debounce the check to avoid excessive calls
    const timer = setTimeout(async () => {
      const startTime = performance.now();
      try {
        const result = await checkPlanarity(nodes, edges);
        const endTime = performance.now();
        const elapsedTime = endTime - startTime;
        setPlanarityCheckTime(elapsedTime);
        setNetworkxPlanarity({ 
          isPlanar: result.isPlanar, 
          checked: result.checked, 
          checking: false,
          graphKey 
        });
      } catch (error) {
        console.error('NetworkX planarity check failed:', error);
        setNetworkxPlanarity({ isPlanar: null, checked: false, checking: false });
        setPlanarityCheckTime(null);
      }
    }, 300); // 300ms debounce

    return () => clearTimeout(timer);
  }, [nodes.length, edges.length, edges]);

  // --------------------------------------------------------------------------
  // DRAG AND DROP HANDLERS - For moving nodes in the interactive graph
  // --------------------------------------------------------------------------

  /**
   * Converts screen coordinates to SVG viewBox coordinates.
   * Uses SVG's built-in matrix transformation for accuracy.
   * Handles both mouse and touch events.
   */
  const getMousePosition = useCallback((e) => {
    const svg = svgRef.current;
    if (!svg) return { x: 0, y: 0 };

    // Handle both mouse and touch events
    const clientX = e.clientX ?? (e.touches?.[0]?.clientX ?? 0);
    const clientY = e.clientY ?? (e.touches?.[0]?.clientY ?? 0);

    // Create a point in screen coordinates
    const point = svg.createSVGPoint();
    point.x = clientX;
    point.y = clientY;

    // Transform to SVG coordinate space
    const ctm = svg.getScreenCTM();
    if (!ctm) return { x: 0, y: 0 };

    const svgPoint = point.matrixTransform(ctm.inverse());
    return { x: svgPoint.x, y: svgPoint.y };
  }, []);

  /**
   * Starts dragging a node when mouse is pressed on it.
   */
  const handleMouseDown = useCallback((e, nodeId) => {
    e.preventDefault();
    e.stopPropagation();
    
    // Note: Branch nodes are prevented from dragging at the JSX level
    // (pointer-events: none and conditional onMouseDown handler)
    // This function will only be called for draggable nodes
    
    // Store initial mouse position and node for drag threshold check
    const pos = getMousePosition(e);
    dragStartPosRef.current = { x: pos.x, y: pos.y };
    draggingNodeRef.current = nodeId;
    isDraggingRef.current = false; // Not dragging yet, just clicked
    setIsDraggingNode(false); // Not dragging yet
    
    // Set offset to 0 so node center stays exactly at mouse position
    dragOffsetRef.current = { x: 0, y: 0 };
  }, [getMousePosition, nodes]);

  /**
   * Updates the dragged node's position as mouse moves.
   * Also handles panning when dragging the canvas.
   * While dragging, node positions are continuously snapped to the nearest
   * available grid point so nodes can never be between grid dots.
   */
  const handleMouseMove = useCallback((e) => {
    // Handle node dragging
    if (draggingNodeRef.current !== null) {
      const pos = getMousePosition(e);
      
      // Check if we've moved enough to start dragging (threshold: 3 pixels)
      if (!isDraggingRef.current) {
        const dx = pos.x - dragStartPosRef.current.x;
        const dy = pos.y - dragStartPosRef.current.y;
        const distance = Math.sqrt(dx * dx + dy * dy);
        
        // Only start dragging if mouse moved more than 3 pixels
        if (distance < 3) {
          return; // Don't move node yet, just a click
        }
        
        // Start dragging now
        isDraggingRef.current = true;
        setIsDraggingNode(true);
      }

      // While dragging, keep snapping node to nearest available grid point
      setNodes(prev => {
        return prev.map(node => {
          if (node.id === draggingNodeRef.current) {
            const snappedPos = snapToNearestAvailable(pos.x, pos.y, prev, node.id);
            return { ...node, x: snappedPos.x, y: snappedPos.y };
          }
          return node;
        });
      });
      return;
    }

    // Handle canvas panning
    if (isPanning2D) {
      const dx = e.clientX - panStart.x;
      const dy = e.clientY - panStart.y;

      // Scale movement based on zoom level
      const scaleFactor = 1 / zoom;
      setPanOffset(prev => ({
        x: prev.x - dx * scaleFactor,
        y: prev.y - dy * scaleFactor
      }));
      setPanStart({ x: e.clientX, y: e.clientY });
    }
  }, [getMousePosition, isPanning2D, panStart, zoom]);

  /**
   * Starts panning when mouse is pressed on canvas background.
   */
  const handleCanvasMouseDown = useCallback((e) => {
    // Only start panning if clicking on the SVG background (not on a node)
    if (e.target.tagName === 'svg' || e.target.classList.contains('grid-background')) {
      setIsPanning2D(true);
      setPanStart({ x: e.clientX, y: e.clientY });
      e.preventDefault();
    }
  }, []);

  /**
   * Ends dragging and recalculates crossings.
   */
  const handleMouseUp = useCallback(() => {
    if (draggingNodeRef.current !== null) {
      // Only snap to grid if we were actually dragging (not just a click)
      if (isDraggingRef.current) {
        // Snap to grid when releasing
        setNodes(currentNodes => {
          const snapped = currentNodes.map(node => {
            if (node.id === draggingNodeRef.current) {
              const snappedPos = snapToNearestAvailable(node.x, node.y, currentNodes, node.id);
              return { ...node, x: snappedPos.x, y: snappedPos.y };
            }
            return node;
          });
          analyzeGraph(snapped, edgesRef.current);
          return snapped;
        });
      }
      
      draggingNodeRef.current = null;
      dragOffsetRef.current = { x: 0, y: 0 };
      isDraggingRef.current = false;
      setIsDraggingNode(false);
    }

    // End panning
    setIsPanning2D(false);
  }, [analyzeGraph]);

  // --------------------------------------------------------------------------
  // ZOOM HANDLERS - For pinch-to-zoom on the graph
  // --------------------------------------------------------------------------

  const lastTouchDistance = useRef(null);

  /**
   * Calculate distance between two touch points.
   */
  const getTouchDistance = (touches) => {
    if (touches.length < 2) return null;
    const dx = touches[0].clientX - touches[1].clientX;
    const dy = touches[0].clientY - touches[1].clientY;
    return Math.sqrt(dx * dx + dy * dy);
  };

  /**
   * Handle touch start for pinch zoom.
   */
  const handleTouchStart = useCallback((e) => {
    if (e.touches.length === 2) {
      e.preventDefault();
      lastTouchDistance.current = getTouchDistance(e.touches);
    }
  }, []);

  /**
   * Handle touch move for pinch zoom.
   */
  const handleTouchMove = useCallback((e) => {
    if (e.touches.length === 2 && lastTouchDistance.current !== null) {
      e.preventDefault();
      const currentDistance = getTouchDistance(e.touches);
      if (currentDistance) {
        const scale = currentDistance / lastTouchDistance.current;
        setZoom(prev => Math.max(0.5, Math.min(3, prev * scale)));
        lastTouchDistance.current = currentDistance;
      }
    }
  }, []);

  /**
   * Handle touch end for pinch zoom.
   */
  const handleTouchEnd = useCallback(() => {
    lastTouchDistance.current = null;
  }, []);

  /**
   * Zoom in button handler.
   */
  const zoomIn = () => setZoom(prev => Math.min(3, prev + 0.2));

  /**
   * Zoom out button handler.
   */
  const zoomOut = () => setZoom(prev => Math.max(0.5, prev - 0.2));

  /**
   * Reset zoom to the default 80% level (zoom factor 0.8).
   */
  const resetZoom = () => {
    setZoom(0.8);
    setPanOffset({ x: 0, y: 0 });
  };

  /**
   * Handle mouse wheel for zooming.
   */
  const handleWheel = useCallback((e) => {
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      const delta = e.deltaY > 0 ? -0.1 : 0.1;
      setZoom(prev => Math.max(0.5, Math.min(3, prev + delta)));
    }
  }, []);

  /**
   * Clear the entire graph.
   */
  const clearGraph = useCallback(() => {
    setNodes([]);
    setEdges([]);
    setCurrentChallenge(null);
    setLoadedExampleId(null);
    setActiveAnalysis(null);  // Close example graph analysis panel
    setShowSolution(false);
    setHighlightedSubgraph(null);
    setShowFixes(false);
    setUploadedFileName(null);
    setMetrics({ total: 0, k: 0, edgeScores: [] });
    setNetworkxPlanarity({ isPlanar: null, checked: false, checking: false });
    resetZoom();
  }, [resetZoom]);

  /**
   * Export graph as JSON.
   */
  const exportGraph = useCallback(() => {
    const graphData = {
      nodes: nodes.map(n => ({ id: n.id, x: n.x, y: n.y })),
      edges: edges.map(e => ({ source: e.source, target: e.target }))
    };
    const blob = new Blob([JSON.stringify(graphData, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `graph-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  }, [nodes, edges]);

  /**
   * Export graph as PNG image.
   */
  const exportAsImage = useCallback(() => {
    if (!svgRef.current || nodes.length === 0) return;
    
    // Create a copy of the SVG with fixed viewBox
    const svg = svgRef.current.cloneNode(true);
    svg.setAttribute('viewBox', '0 0 400 400');
    svg.setAttribute('width', '800');
    svg.setAttribute('height', '800');
    
    const svgData = new XMLSerializer().serializeToString(svg);
    const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(svgBlob);
    
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = 800;
      canvas.height = 800;
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#0a0f1c';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      
      canvas.toBlob((blob) => {
        if (blob) {
          const downloadUrl = URL.createObjectURL(blob);
          const a = document.createElement('a');
          a.href = downloadUrl;
          a.download = `graph-${Date.now()}.png`;
          document.body.appendChild(a);
          a.click();
          document.body.removeChild(a);
          URL.revokeObjectURL(downloadUrl);
        }
        URL.revokeObjectURL(url);
      }, 'image/png');
    };
    img.onerror = () => {
      console.error('Failed to export image');
      URL.revokeObjectURL(url);
    };
    img.src = url;
  }, [nodes.length]);

  // Keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (e) => {
      // Ctrl/Cmd + Z: Clear graph
      if ((e.ctrlKey || e.metaKey) && e.key === 'z' && !e.shiftKey) {
        e.preventDefault();
        if (e.ctrlKey || e.metaKey) {
          clearGraph();
        }
      }
      // Delete/Backspace: Clear graph when focused
      if ((e.key === 'Delete' || e.key === 'Backspace') && !e.target.matches('input, textarea')) {
        if (selectedNode) {
          // Remove selected node and its edges
          setNodes(prev => prev.filter(n => n.id !== selectedNode));
          setEdges(prev => prev.filter(e => e.source !== selectedNode && e.target !== selectedNode));
          setSelectedNode(null);
        } else if (nodes.length > 0) {
          clearGraph();
        }
      }
      // Escape: Deselect
      if (e.key === 'Escape') {
        setSelectedNode(null);
        setHoveredEdge(null);
      }
      // G: Toggle grid
      if (e.key === 'g' && !e.target.matches('input, textarea')) {
        e.preventDefault();
        setShowGrid(prev => {
          if (prev) {
            return false; // Turn off
          } else {
            // Turn on with last used style (or 'dots' as default)
            return true;
          }
        });
      }
      // 0: Reset zoom
      if (e.key === '0' && !e.target.matches('input, textarea')) {
        e.preventDefault();
        resetZoom();
      }
    };

    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [clearGraph, resetZoom, selectedNode, nodes.length]);

  /**
   * Calculate the viewBox based on zoom level.
   */
  const getViewBox = () => {
    const size = 400 / zoom;
    const offset = (400 - size) / 2;
    return `${offset + panOffset.x} ${offset + panOffset.y} ${size} ${size}`;
  };

  // --------------------------------------------------------------------------
  // FILE UPLOAD HANDLERS - For importing graph files
  // --------------------------------------------------------------------------

  /**
   * Loads a parsed graph into the interactive playground.
   * Clears any active challenge and resets the solution display.
   */
  const loadGraph = useCallback((graphNodes, graphEdges, fileName) => {
    setNodes(graphNodes);
    setEdges(graphEdges);
    setCurrentChallenge(null);
    setLoadedExampleId(null);  // Clear any loaded example
    setActiveAnalysis(null);  // Close example graph analysis panel
    setShowSolution(false);
    setHighlightedSubgraph(null);
    setShowFixes(false);
    setUploadedFileName(fileName);
    analyzeGraph(graphNodes, graphEdges);
    setUploadStatus({
      type: 'success',
      message: `Loaded ${graphNodes.length} nodes and ${graphEdges.length} edges`
    });

    // Clear status after 5 seconds
    setTimeout(() => setUploadStatus(null), 5000);
  }, [analyzeGraph]);

  /**
   * Processes an uploaded file - reads and parses it.
   */
  const processFile = useCallback(async (file) => {
    // If an upload is already in progress, ignore additional requests to avoid
    // stacking heavy work on the main thread.
    if (isUploading) return;

    setIsUploading(true);

    // Validate file type
    const validTypes = ['application/json', 'text/xml', 'application/xml', ''];
    const validExtensions = ['json', 'graphml', 'xml'];
    const extension = file.name.split('.').pop().toLowerCase();

    if (!validTypes.includes(file.type) && !validExtensions.includes(extension)) {
      setUploadStatus({
        type: 'error',
        message: 'Invalid file type. Please upload .json or .graphml files.'
      });
      return;
    }

    try {
      const content = await file.text();
      const { nodes: parsedNodes, edges: parsedEdges } = parseGraphFile(file, content);

      if (!Array.isArray(parsedNodes) || parsedNodes.length === 0) {
        throw new Error('No nodes found in the file.');
      }

      if (!Array.isArray(parsedEdges)) {
        throw new Error('File format error: edges array is missing or invalid.');
      }

      // Hard reject extremely large graphs before we try to render or analyze
      // them. This avoids React rendering and layout work on thousands of
      // elements which can lock up the browser.
      if (parsedNodes.length > MAX_UPLOAD_NODES || parsedEdges.length > MAX_UPLOAD_EDGES) {
        setUploadStatus({
          type: 'error',
          message: `Uploaded graph is too large to display safely (nodes: ${parsedNodes.length}, edges: ${parsedEdges.length}).` +
                   ` Please upload a smaller graph or simplify it before importing.`
        });
        // Auto-clear the message after a few seconds
        setTimeout(() => setUploadStatus(null), 6000);
        return;
      }

      loadGraph(parsedNodes, parsedEdges, file.name);
    } catch (error) {
      console.error('File parsing error:', error);
      setUploadStatus({
        type: 'error',
        message: error.message || 'Failed to parse file. Please check the format.'
      });

      // Clear error after 5 seconds
      setTimeout(() => setUploadStatus(null), 5000);
    } finally {
      setIsUploading(false);
    }
  }, [isUploading, loadGraph]);

  const handleDragOver = (e) => {
    // Prevent the browser from treating this as a navigation and
    // indicate that we're actively handling the drop.
    e.preventDefault();
    e.stopPropagation();
    setIsDragOver(true);
  };

  const handleDragLeave = () => {
    setIsDragOver(false);
  };

  const handleDrop = (e) => {
    // Prevent default browser behaviour which would try to open the
    // dropped file in the tab (making the app appear to "crash").
    e.preventDefault();
    e.stopPropagation();
    setIsDragOver(false);

    const file = e.dataTransfer.files?.[0];
    if (file) {
      processFile(file);
    }
  };

  // Prevent full-page navigation when files are dragged/dropped anywhere
  // on the window instead of just inside our dedicated drop zone.
  useEffect(() => {
    const preventWindowFileDrop = (e) => {
      e.preventDefault();
      e.stopPropagation();
    };

    window.addEventListener('dragover', preventWindowFileDrop);
    window.addEventListener('drop', preventWindowFileDrop);

    return () => {
      window.removeEventListener('dragover', preventWindowFileDrop);
      window.removeEventListener('drop', preventWindowFileDrop);
    };
  }, []);

  /**
   * Triggers the hidden file input when "Browse Files" button is clicked.
   */
  const handleBrowseClick = () => {
    fileInputRef.current?.click();
  };

  /**
   * Handles file selection from the file input.
   */
  const handleFileSelect = (e) => {
    const file = e.target.files?.[0];
    if (file) {
      processFile(file);
    }
    // Reset input so the same file can be selected again
    e.target.value = '';
  };

  // --------------------------------------------------------------------------
  // EXAMPLE GRAPH ANALYSIS - For the "Analyze" button on example cards
  // --------------------------------------------------------------------------

  /**
   * Generates and displays the analysis panel for an example graph.
   */
  const handleAnalyzeGraph = (graphId) => {
    const graphData = EXAMPLE_GRAPHS[graphId];
    if (graphData) {
      const { nodes: graphNodes, edges: graphEdges } = graphData.generate();
      const crossings = calculateCrossings(graphNodes, graphEdges);
      setActiveAnalysis({
        ...graphData,
        nodes: graphNodes,
        edges: graphEdges,
        crossings: crossings.total,
        edgeScores: crossings.edgeScores
      });
    }
  };

  /**
   * Loads an example graph onto the grid (without challenge mode).
   */
  const loadExampleGraph = (graphId) => {
    const graphData = EXAMPLE_GRAPHS[graphId];
    if (graphData) {
      const { nodes: graphNodes, edges: graphEdges } = graphData.generate();
      setNodes(graphNodes);
      setEdges(graphEdges);
      setCurrentChallenge(null);  // Not a challenge, just loading
      setLoadedExampleId(graphId);  // Track which example was loaded
      setShowSolution(false);
      setHighlightedSubgraph(null);
      setShowFixes(false);
      setOriginalPositions(null);
      analyzeGraph(graphNodes, graphEdges);
    }
  };

  /**
   * Closes the analysis panel.
   */
  const closeAnalysis = () => {
    setActiveAnalysis(null);
  };

  // --------------------------------------------------------------------------
  // CHALLENGE MODE - For loading example graphs into the interactive playground
  // --------------------------------------------------------------------------

  /**
   * Loads an example graph into the interactive playground as a challenge.
   */
  const loadChallengeGraph = (graphId) => {
    const graphData = EXAMPLE_GRAPHS[graphId];
    if (graphData) {
      const { nodes: graphNodes, edges: graphEdges } = graphData.generate();
      setNodes(graphNodes);
      setEdges(graphEdges);
      setCurrentChallenge(graphId);
      setShowSolution(false);  // Hide solution when loading new challenge
      setHighlightedSubgraph(null);  // Clear any highlighting
      setShowFixes(false);
      analyzeGraph(graphNodes, graphEdges);
    }
  };

  // --------------------------------------------------------------------------
  // SOLUTION LOGIC - Determines if a graph is solvable (planar)
  // --------------------------------------------------------------------------

  /**
   * Determines whether the current graph is solvable (planar) and provides explanation.
   * Uses Euler's formula (E ≤ 3V - 6) as a necessary condition for planarity.
   * Also detects K5 and K3,3 subgraphs.
   */
  const getSolutionInfo = useCallback(() => {
    const thoughtProcess = [];
    
    const v = nodes.length;  // Number of vertices
    const e = edges.length;  // Number of edges
    
    // Step 1: Check graph size
    thoughtProcess.push({
      step: 1,
      description: `Check graph size: ${v} vertices, ${e} edges`,
      result: v < 3 ? 'Graph has fewer than 3 vertices - always planar' : 'Graph has 3+ vertices - proceed with analysis'
    });
    
    if (v < 3) {
      return { 
        solvable: true, 
        reason: 'Graphs with fewer than 3 vertices are always planar.', 
        k5: null, 
        k33: null, 
        famousGraph: null, 
        fixes: null,
        thoughtProcess
      };
    }

    // For very large graphs we STILL attempt to find K₅ / K₃,₃ witnesses,
    // but we keep the most expensive "minimum removal" searches guarded
    // elsewhere. This respects Kuratowski's theorem while keeping the UI safe.
    if (v > MAX_NODES_FOR_MIN_REMOVAL || e > MAX_EDGES_FOR_MIN_REMOVAL) {
      thoughtProcess.push({
        step: 2,
        description: `Graph is large: ${v} vertices, ${e} edges`,
        result: 'Skipping minimum-removal search, but still hunting for K₅ / K₃,₃ subgraphs.'
      });
    }
    
    // Check for famous graph first
    const famousGraphId = identifyFamousGraph(nodes, edges);
    let famousGraphInfo = null;
    if (famousGraphId && EXAMPLE_GRAPHS[famousGraphId]) {
      famousGraphInfo = EXAMPLE_GRAPHS[famousGraphId];
      thoughtProcess.push({
        step: 2,
        description: `Check for famous graph identification`,
        result: `Identified as: ${famousGraphInfo.name}`
      });
    } else {
      thoughtProcess.push({
        step: 2,
        description: `Check for famous graph identification`,
        result: 'Not a known famous graph'
      });
    }

    // Check if this is a known non-planar example graph (either as challenge or loaded)
    const knownGraphId = currentChallenge || loadedExampleId;
    if (knownGraphId && EXAMPLE_GRAPHS[knownGraphId]) {
      const graphData = EXAMPLE_GRAPHS[knownGraphId];
      thoughtProcess.push({
        step: 3,
        description: `Check if graph matches known example: ${knownGraphId}`,
        result: `Matched known graph: ${graphData.name} (${graphData.isPlanar ? 'planar' : 'non-planar'})`
      });

      // For K5, highlight all 5 nodes
      let k5Info = null;
      if (knownGraphId === 'k5') {
        k5Info = { nodes: nodes.map(n => n.id), edges };
      }

      // For K3,3, highlight with proper partition
      let k33Info = null;
      if (knownGraphId === 'k33') {
        k33Info = { nodes: nodes.map(n => n.id), edges, setA: [0, 1, 2], setB: [3, 4, 5] };
      }

      // For graphs that contain K5 or K3,3 as minors/subdivisions, try to find them
      // Petersen graph contains K5 as a minor (but not as subgraph)
      // Heawood, Pappus, etc. contain K3,3 as subdivision
      if (!k5Info && !k33Info) {
        // Try exact subgraph detection first
        k5Info = findK5Subgraph(nodes, edges);
        k33Info = findK33Subgraph(nodes, edges);
        
        // If no exact subgraphs found, try subdivision detection
        if (!k5Info && !k33Info) {
          k5Info = findK5Subdivision(nodes, edges);
          k33Info = findK33Subdivision(nodes, edges);
        }
      }

      // Find fixes if non-planar
      let fixes = null;
      if (graphData.isPlanar === false) {
        fixes = findPlanarityFixes(nodes, edges);
      }

      thoughtProcess.push({
        step: 4,
        description: 'Final conclusion',
        result: graphData.isPlanar !== false ? 'Graph is planar' : 'Graph is non-planar (known example)'
      });

      return {
        solvable: graphData.isPlanar !== false ? true : false,
        reason: graphData.explanation,
        k5: k5Info,
        k33: k33Info,
        famousGraph: famousGraphInfo,
        fixes: fixes,
        thoughtProcess
      };
    }

    // Step 3: Check NetworkX result
    thoughtProcess.push({
      step: 3,
      description: 'Check NetworkX planarity result',
      result: networkxPlanarity.checked 
        ? (networkxPlanarity.isPlanar !== null 
          ? `NetworkX result: ${networkxPlanarity.isPlanar ? 'planar' : 'non-planar'}` 
          : 'NetworkX check completed but result is null')
        : (networkxPlanarity.checking ? 'NetworkX is currently checking...' : 'NetworkX not yet checked')
    });

    // Try to find K5 or K3,3 subgraphs for display purposes
    let k5 = findK5Subgraph(nodes, edges);
    let k33 = findK33Subgraph(nodes, edges);
    
    // If no exact subgraphs found, try subdivision detection
    if (!k5 && !k33) {
      k5 = findK5Subdivision(nodes, edges);
      k33 = findK33Subdivision(nodes, edges);
    }

    // Step 4: Check for K5 and K3,3
    thoughtProcess.push({
      step: 4,
      description: 'Check for K₅ and K₃,₃ subgraphs/subdivisions',
      result: k5 && k33 
        ? `Found both K₅ (${k5.branchNodes ? 'subdivision' : 'subgraph'}) and K₃,₃ (${k33.branchNodes ? 'subdivision' : 'subgraph'})`
        : k5 
          ? `Found K₅ as ${k5.branchNodes ? 'subdivision' : 'subgraph'}`
          : k33 
            ? `Found K₃,₃ as ${k33.branchNodes ? 'subdivision' : 'subgraph'}`
            : 'No K₅ or K₃,₃ found'
    });

    // Use NetworkX result if available (definitive answer)
    if (networkxPlanarity.checked && networkxPlanarity.isPlanar !== null) {
      if (networkxPlanarity.isPlanar) {
        // Graph IS planar according to NetworkX
        thoughtProcess.push({
          step: 5,
          description: `Check crossing count: ${metrics.total} crossings`,
          result: metrics.total === 0 ? 'No crossings - perfect planar embedding' : `${metrics.total} crossings present but graph is planar`
        });
        thoughtProcess.push({
          step: 6,
          description: 'Final conclusion',
          result: 'Graph is PLANAR (verified by NetworkX)'
        });
        
        if (metrics.total === 0) {
          return { 
            solvable: true, 
            reason: 'This graph is planar! You\'ve found a valid planar embedding with no crossings.', 
            k5: null, 
            k33: null, 
            famousGraph: famousGraphInfo, 
            fixes: null,
            thoughtProcess
          };
        }
        return { 
          solvable: true, 
          reason: `This graph IS planar (verified by NetworkX). It has ${v} vertices and ${e} edges. A planar embedding exists — try the Untangle button or manually eliminate the ${metrics.total} crossing${metrics.total > 1 ? 's' : ''}!`, 
          k5: null, 
          k33: null, 
          famousGraph: famousGraphInfo, 
          fixes: null,
          thoughtProcess
        };
      } else {
        // Graph is NOT planar according to NetworkX
        let reason = `This graph is NOT planar (verified by NetworkX). It has ${v} vertices and ${e} edges.`;
        if (k5 && k33) {
          const k5Type = k5.branchNodes ? 'subdivision' : 'subgraph';
          const k33Type = k33.branchNodes ? 'subdivision' : 'subgraph';
          reason += ` Contains both K₅ (as ${k5Type}) and K₃,₃ (as ${k33Type}).`;
        } else if (k5) {
          const k5Type = k5.branchNodes ? 'subdivision' : 'subgraph';
          reason += ` Contains K₅ as a ${k5Type}.`;
        } else if (k33) {
          const k33Type = k33.branchNodes ? 'subdivision' : 'subgraph';
          reason += ` Contains K₃,₃ as a ${k33Type}.`;
        } else {
          reason += ' It contains a K₅ or K₃,₃ subdivision (Kuratowski\'s theorem).';
        }

        if (famousGraphInfo) {
          reason = `Identified as ${famousGraphInfo.name}. ${famousGraphInfo.explanation}`;
        }

        thoughtProcess.push({
          step: 5,
          description: 'Final conclusion',
          result: 'Graph is NON-PLANAR (verified by NetworkX)'
        });

        const fixes = findPlanarityFixes(nodes, edges);
        return { solvable: false, reason, k5, k33, famousGraph: famousGraphInfo, fixes, thoughtProcess };
      }
    }

    // NetworkX is still checking or hasn't loaded yet - show loading state
    if (networkxPlanarity.checking) {
      thoughtProcess.push({
        step: 5,
        description: 'Final conclusion',
        result: 'Waiting for NetworkX to complete check...'
      });
      return { 
        solvable: true, 
        reason: `Checking planarity with NetworkX... (${v} vertices, ${e} edges)`, 
        k5: null, 
        k33: null, 
        famousGraph: famousGraphInfo, 
        fixes: null,
        loading: true,
        thoughtProcess
      };
    }

    // Fallback: NetworkX not available, use heuristics
    // Euler's formula check: E ≤ 3V - 6 is necessary for planarity
    thoughtProcess.push({
      step: 5,
      description: `Check Euler's formula: E ≤ 3V - 6 (${e} ≤ ${3 * v - 6})`,
      result: e > 3 * v - 6 
        ? `Violates Euler's formula: ${e} > ${3 * v - 6} - cannot be planar`
        : `Satisfies Euler's formula: ${e} ≤ ${3 * v - 6} - proceed with further checks`
    });

    if (e > 3 * v - 6) {
      let reason = `This graph has ${e} edges but a planar graph with ${v} vertices can have at most ${3 * v - 6} edges (Euler's formula). It cannot be planar.`;
      if (k5) {
        const k5Type = k5.branchNodes ? 'subdivision' : 'subgraph';
        reason += ` Contains K₅ as a ${k5Type}.`;
      }
      if (k33) {
        const k33Type = k33.branchNodes ? 'subdivision' : 'subgraph';
        reason += ` Contains K₃,₃ as a ${k33Type}.`;
      }

      thoughtProcess.push({
        step: 6,
        description: 'Final conclusion',
        result: 'Graph is NON-PLANAR (violates Euler\'s formula)'
      });

      const fixes = findPlanarityFixes(nodes, edges);
      return { solvable: false, reason, k5, k33, famousGraph: famousGraphInfo, fixes, thoughtProcess };
    }

    // Check if K5 or K3,3 was found
    if (k5 || k33) {
      let reason = 'This graph contains ';
      if (k5 && k33) {
        const k5Type = k5.branchNodes ? 'subdivision' : 'subgraph';
        const k33Type = k33.branchNodes ? 'subdivision' : 'subgraph';
        reason += `both K₅ (as ${k5Type}) and K₃,₃ (as ${k33Type})`;
      } else if (k5) {
        const k5Type = k5.branchNodes ? 'subdivision' : 'subgraph';
        reason += `K₅ (complete graph on 5 vertices) as a ${k5Type}`;
      } else {
        const k33Type = k33.branchNodes ? 'subdivision' : 'subgraph';
        reason += `K₃,₃ (complete bipartite graph) as a ${k33Type}`;
      }
      reason += '. By Kuratowski\'s theorem, it cannot be planar.';

      thoughtProcess.push({
        step: 6,
        description: 'Final conclusion',
        result: 'Graph is NON-PLANAR (contains K₅ or K₃,₃)'
      });

      const fixes = findPlanarityFixes(nodes, edges);
      return { solvable: false, reason, k5, k33, famousGraph: famousGraphInfo, fixes, thoughtProcess };
    }

    // If already at 0 crossings, it's definitely planar
    thoughtProcess.push({
      step: 5,
      description: `Check crossing count: ${metrics.total} crossings`,
      result: metrics.total === 0 ? 'No crossings - graph is planar' : `${metrics.total} crossings present`
    });

    if (metrics.total === 0) {
      thoughtProcess.push({
        step: 6,
        description: 'Final conclusion',
        result: 'Graph is PLANAR (no crossings in current embedding)'
      });
      return { solvable: true, reason: 'This graph is planar! You\'ve already found an embedding with no crossings.', k5: null, k33: null, famousGraph: famousGraphInfo, fixes: null, thoughtProcess };
    }

    // Heuristics inconclusive, waiting for NetworkX
    thoughtProcess.push({
      step: 6,
      description: 'Final conclusion',
      result: 'Heuristics inconclusive - waiting for NetworkX check'
    });
    return {
      solvable: true,
      reason: `Analyzing planarity... (${v} vertices, ${e} edges, ${metrics.total} crossings). Click Untangle to check with NetworkX.`,
      k5: null,
      k33: null,
      famousGraph: famousGraphInfo,
      fixes: null,
      thoughtProcess
    };
  }, [nodes, edges, metrics.total, currentChallenge, loadedExampleId, networkxPlanarity]);

  /**
   * Gets the edges that should be highlighted based on the selected subgraph.
   * For subdivisions, only highlights edges that are part of the subdivision paths.
   */
  const getHighlightedEdges = useCallback(() => {
    if (showFixes) {
      const solutionInfo = getSolutionInfo();
      if (solutionInfo.fixes && solutionInfo.fixes.edges) {
        return solutionInfo.fixes.edges;
      }
      return new Set();
    }
    if (!highlightedSubgraph) return new Set();

    const solutionInfo = getSolutionInfo();
    const edgeSet = new Set();

    if (highlightedSubgraph === 'k5' && solutionInfo.k5) {
      if (solutionInfo.k5.branchNodes && solutionInfo.k5.paths) {
        // It's a subdivision - only highlight edges in the paths
        // paths is a Map, so iterate over values
        for (const path of solutionInfo.k5.paths.values()) {
          for (let i = 0; i < path.length - 1; i++) {
            const source = path[i];
            const target = path[i + 1];
            edgeSet.add(`${Math.min(source, target)}-${Math.max(source, target)}`);
          }
        }
      } else {
        // It's an exact subgraph - highlight all edges
        (solutionInfo.k5.edges || []).forEach(e => {
          edgeSet.add(`${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`);
        });
      }
    } else if (highlightedSubgraph === 'k33' && solutionInfo.k33) {
      if (solutionInfo.k33.branchNodes && solutionInfo.k33.paths) {
        // It's a subdivision - only highlight edges in the paths
        // paths is a Map, so iterate over values
        for (const path of solutionInfo.k33.paths.values()) {
          for (let i = 0; i < path.length - 1; i++) {
            const source = path[i];
            const target = path[i + 1];
            edgeSet.add(`${Math.min(source, target)}-${Math.max(source, target)}`);
          }
        }
      } else {
        // It's an exact subgraph - highlight all edges
        (solutionInfo.k33.edges || []).forEach(e => {
          edgeSet.add(`${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`);
        });
      }
    }

    return edgeSet;
  }, [highlightedSubgraph, showFixes, getSolutionInfo]);

  /**
   * Gets the nodes that should be highlighted based on the selected subgraph.
   * For subdivisions, only returns branch nodes (5 for K5, 6 for K3,3).
   * For exact subgraphs, returns all nodes in the subgraph.
   */
  const getHighlightedNodes = useCallback(() => {
    if (showFixes) {
      const solutionInfo = getSolutionInfo();
      if (solutionInfo.fixes && solutionInfo.fixes.nodes) {
        return solutionInfo.fixes.nodes;
      }
      return new Set();
    }
    if (!highlightedSubgraph) return new Set();

    const solutionInfo = getSolutionInfo();

    if (highlightedSubgraph === 'k5' && solutionInfo.k5) {
      // For subdivisions, only highlight branch nodes (5 nodes)
      if (solutionInfo.k5.branchNodes) {
        return new Set(solutionInfo.k5.branchNodes);
      }
      // For exact subgraphs, return all nodes
      return new Set(solutionInfo.k5.nodes || []);
    } else if (highlightedSubgraph === 'k33' && solutionInfo.k33) {
      // For subdivisions, only highlight branch nodes (6 nodes)
      if (solutionInfo.k33.branchNodes) {
        const branchNodes = [...(solutionInfo.k33.branchNodes.setA || []), ...(solutionInfo.k33.branchNodes.setB || [])];
        return new Set(branchNodes);
      }
      // For exact subgraphs, return all nodes
      return new Set(solutionInfo.k33.nodes || []);
    }

    return new Set();
  }, [highlightedSubgraph, showFixes, getSolutionInfo]);

  /**
   * Gets the branch nodes (main vertices) of a subdivision for special highlighting.
   * This is the same as getHighlightedNodes for subdivisions, but kept for clarity.
   */
  const getBranchNodes = useCallback(() => {
    return getHighlightedNodes();
  }, [getHighlightedNodes]);

  /**
   * Gets subdivision nodes (intermediate nodes on paths) for less prominent highlighting.
   */
  const getSubdivisionNodes = useCallback(() => {
    if (!highlightedSubgraph || showFixes) return new Set();
    
    const solutionInfo = getSolutionInfo();
    const subdivisionNodes = new Set();
    const branchNodeSet = getHighlightedNodes();
    
    if (highlightedSubgraph === 'k5' && solutionInfo.k5 && solutionInfo.k5.branchNodes && solutionInfo.k5.paths) {
      // Collect all nodes in paths that are not branch nodes
      for (const path of solutionInfo.k5.paths.values()) {
        for (const nodeId of path) {
          if (!branchNodeSet.has(nodeId)) {
            subdivisionNodes.add(nodeId);
          }
        }
      }
    } else if (highlightedSubgraph === 'k33' && solutionInfo.k33 && solutionInfo.k33.branchNodes && solutionInfo.k33.paths) {
      // Collect all nodes in paths that are not branch nodes
      for (const path of solutionInfo.k33.paths.values()) {
        for (const nodeId of path) {
          if (!branchNodeSet.has(nodeId)) {
            subdivisionNodes.add(nodeId);
          }
        }
      }
    }
    
    return subdivisionNodes;
  }, [highlightedSubgraph, showFixes, getSolutionInfo, getHighlightedNodes]);

  /**
   * Calculates 3D sphere positions for nodes using Fibonacci sphere distribution.
   * Projects 3D coordinates to 2D for rendering with rotation support.
   */
  const get3DNodePositions = useCallback(() => {
    if (nodes.length === 0) return [];

    const n = nodes.length;
    const radius = 120;
    const cx = 200, cy = 200;
    const goldenRatio = (1 + Math.sqrt(5)) / 2;

    // Rotation angles in radians
    const rotX = sphereRotation.x * Math.PI / 180;
    const rotY = sphereRotation.y * Math.PI / 180;

    return nodes.map((node, i) => {
      // Fibonacci sphere point distribution
      const theta = 2 * Math.PI * i / goldenRatio;
      const phi = Math.acos(1 - 2 * (i + 0.5) / n);

      // 3D coordinates on unit sphere
      let x = Math.sin(phi) * Math.cos(theta);
      let y = Math.sin(phi) * Math.sin(theta);
      let z = Math.cos(phi);

      // Apply Y rotation (left-right)
      const cosY = Math.cos(rotY);
      const sinY = Math.sin(rotY);
      const x2 = x * cosY + z * sinY;
      const z2 = -x * sinY + z * cosY;
      x = x2;
      z = z2;

      // Apply X rotation (up-down)
      const cosX = Math.cos(rotX);
      const sinX = Math.sin(rotX);
      const y2 = y * cosX - z * sinX;
      const z3 = y * sinX + z * cosX;
      y = y2;
      z = z3;

      // Perspective projection (closer = larger)
      const perspective = 400;
      const scale = perspective / (perspective + z * radius);

      return {
        id: node.id,
        x: cx + x * radius * scale,
        y: cy + y * radius * scale,
        z: z,
        scale: scale,
        opacity: 0.4 + 0.6 * ((z + 1) / 2) // Front nodes more opaque
      };
    });
  }, [nodes, sphereRotation]);

  /**
   * Calculates edge positions for 3D view with depth sorting.
   */
  const get3DEdges = useCallback(() => {
    const positions = get3DNodePositions();
    const posMap = new Map(positions.map(p => [p.id, p]));

    return edges.map(edge => {
      const start = posMap.get(edge.source);
      const end = posMap.get(edge.target);
      if (!start || !end) return null;

      const avgZ = (start.z + end.z) / 2;
      return {
        ...edge,
        x1: start.x, y1: start.y,
        x2: end.x, y2: end.y,
        z: avgZ,
        opacity: 0.3 + 0.5 * ((avgZ + 1) / 2)
      };
    }).filter(Boolean).sort((a, b) => a.z - b.z); // Back edges first
  }, [edges, get3DNodePositions]);

  // Auto-rotate the sphere when 3D view is active (pause when dragging or disabled)
  useEffect(() => {
    if (!show3DView || isDragging3D || !autoRotate3D) return;

    const interval = setInterval(() => {
      setSphereRotation(prev => ({
        x: prev.x,
        y: (prev.y + 0.5) % 360
      }));
    }, 50);

    return () => clearInterval(interval);
  }, [show3DView, isDragging3D, autoRotate3D]);


  /**
   * Calculates canonical positions for K5 (pentagon) or K3,3 (bipartite) layout.
   */
  const getCanonicalPositions = useCallback((type, nodeIds, k33Info) => {
    const cx = 200, cy = 200; // Center of view

    if (type === 'k5' && nodeIds.length === 5) {
      // Pentagon layout for K5
      const radius = 100;
      const positions = {};
      for (let i = 0; i < 5; i++) {
        const angle = (i * 2 * Math.PI / 5) - Math.PI / 2; // Start from top
        positions[nodeIds[i]] = {
          x: cx + radius * Math.cos(angle),
          y: cy + radius * Math.sin(angle)
        };
      }
      return positions;
    }

    if (type === 'k33' && k33Info && k33Info.setA && k33Info.setB) {
      // Bipartite layout for K3,3
      const positions = {};
      const leftX = cx - 80;
      const rightX = cx + 80;
      const yStart = cy - 80;
      const yStep = 80;

      k33Info.setA.forEach((nodeId, i) => {
        positions[nodeId] = { x: leftX, y: yStart + i * yStep };
      });
      k33Info.setB.forEach((nodeId, i) => {
        positions[nodeId] = { x: rightX, y: yStart + i * yStep };
      });
      return positions;
    }

    return null;
  }, []);

  /**
   * Animates a single node to its target position smoothly.
   * Uses a ref to prevent multiple simultaneous animations.
   */
  const animateSingleNode = useCallback((nodeId, startPos, targetPos, duration) => {
    return new Promise((resolve) => {
      let animationFrameId = null;
      let cancelled = false;
      
      const startTime = performance.now();
      setCurrentlyMovingNode(nodeId); // Mark this node as currently moving

      const animate = (currentTime) => {
        if (cancelled) {
          resolve();
          return;
        }
        
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);

        // Easing function (ease-in-out cubic for smooth movement)
        const eased = progress < 0.5
          ? 4 * progress * progress * progress
          : 1 - Math.pow(-2 * progress + 2, 3) / 2;

        if (progress < 1) {
          setNodes(prevNodes => prevNodes.map(node => {
            if (node.id === nodeId) {
              const newX = startPos.x + (targetPos.x - startPos.x) * eased;
              const newY = startPos.y + (targetPos.y - startPos.y) * eased;
              return {
                ...node,
                x: newX,
                y: newY
              };
            }
            return node;
          }));
          animationFrameId = requestAnimationFrame(animate);
        } else {
          // Final update to ensure exact position - only once
          setNodes(prevNodes => prevNodes.map(node => {
            if (node.id === nodeId) {
              return {
                ...node,
                x: targetPos.x,
                y: targetPos.y
              };
            }
            return node;
          }));
          setCurrentlyMovingNode(null); // Clear moving indicator
          resolve();
        }
      };

      animationFrameId = requestAnimationFrame(animate);
      
      // Return cleanup function
      return () => {
        cancelled = true;
        if (animationFrameId) {
          cancelAnimationFrame(animationFrameId);
        }
      };
    });
  }, []);

  /**
   * Animates nodes to their target positions one by one using a single RAF loop.
   * This is smoother than async/await because it updates only one node per frame.
   * For subdivisions, only animates branch nodes - subdivision nodes stay in place.
   */
  const animateToPositionsSequential = useCallback((targetPositions, perNodeDuration = 1200, pauseBetween = 300) => {
    if (!targetPositions) return;

    setIsAnimating(true);

    // Get the node IDs that need to move
    const nodeIdsToMove = Object.keys(targetPositions).map(id =>
      typeof id === 'string' ? parseInt(id, 10) : id
    );

    // Capture starting positions from current nodes
    const startPositions = {};
    nodes.forEach(node => {
      if (targetPositions[node.id] !== undefined) {
        startPositions[node.id] = { x: node.x, y: node.y };
      }
    });

    let animationFrameId = null;
    let currentNodeIndex = 0;
    let currentPhase = 'moving'; // 'moving' or 'pausing'
    let phaseStartTime = performance.now();
    let currentNodeId = nodeIdsToMove[currentNodeIndex];
    let currentNodeStart = startPositions[currentNodeId];
    let currentNodeTarget = targetPositions[currentNodeId];

    const animate = (currentTime) => {
      if (currentNodeIndex >= nodeIdsToMove.length) {
        // All nodes animated - final update to ensure exact positions
        setNodes(prevNodes => {
          const updated = prevNodes.map(node => {
            if (targetPositions[node.id] !== undefined) {
              const targetPos = targetPositions[node.id];
              return {
                ...node,
                x: targetPos.x,
                y: targetPos.y
              };
            }
            return node;
          });
          return updated;
        });
        setIsAnimating(false);
        return;
      }

      const elapsed = currentTime - phaseStartTime;

      if (currentPhase === 'moving') {
        // Animate current node
        const progress = Math.min(elapsed / perNodeDuration, 1);
        
        // Simpler easing function for better performance (ease-in-out quadratic)
        const eased = progress < 0.5
          ? 2 * progress * progress
          : 1 - Math.pow(-2 * progress + 2, 2) / 2;

        if (progress < 1) {
          // Update only the current node
          setNodes(prevNodes => prevNodes.map(node => {
            if (node.id === currentNodeId && currentNodeStart && currentNodeTarget) {
              return {
                ...node,
                x: currentNodeStart.x + (currentNodeTarget.x - currentNodeStart.x) * eased,
                y: currentNodeStart.y + (currentNodeTarget.y - currentNodeStart.y) * eased
              };
            }
            return node;
          }));
          animationFrameId = requestAnimationFrame(animate);
        } else {
          // Current node finished - move to next or pause
          // Ensure exact final position for current node
          setNodes(prevNodes => prevNodes.map(node => {
            if (node.id === currentNodeId && currentNodeTarget) {
              return {
                ...node,
                x: currentNodeTarget.x,
                y: currentNodeTarget.y
              };
            }
            return node;
          }));

          // Move to next node
          currentNodeIndex++;
          if (currentNodeIndex < nodeIdsToMove.length) {
            // Start pause phase
            currentPhase = 'pausing';
            phaseStartTime = currentTime;
            animationFrameId = requestAnimationFrame(animate);
          } else {
            // All done
            setIsAnimating(false);
          }
        }
      } else if (currentPhase === 'pausing') {
        // Pause between nodes
        if (elapsed >= pauseBetween) {
          // Pause done - start next node
          currentNodeId = nodeIdsToMove[currentNodeIndex];
          currentNodeStart = startPositions[currentNodeId];
          currentNodeTarget = targetPositions[currentNodeId];
          currentPhase = 'moving';
          phaseStartTime = currentTime;
        }
        animationFrameId = requestAnimationFrame(animate);
      }
    };

    animationFrameId = requestAnimationFrame(animate);

    // Return cleanup function
    return () => {
      if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
      }
      setIsAnimating(false);
    };
  }, [nodes]);

  /**
   * Animates all nodes back to original positions simultaneously (faster for return).
   */
  const animateToPositionsSimultaneous = useCallback((targetPositions, duration = 500) => {
    if (!targetPositions) return;

    setIsAnimating(true);
    let animationFrameId = null;
    const startTime = performance.now();
    const startPositions = {};

    // Capture starting positions - use ref to avoid stale closure
    const currentNodes = nodes;
    currentNodes.forEach(node => {
      if (targetPositions[node.id] !== undefined) {
        startPositions[node.id] = { x: node.x, y: node.y };
      }
    });

    const animate = (currentTime) => {
      const elapsed = currentTime - startTime;
      const progress = Math.min(elapsed / duration, 1);

      // Simpler easing function for better performance (ease-in-out quadratic)
      const eased = progress < 0.5
        ? 2 * progress * progress
        : 1 - Math.pow(-2 * progress + 2, 2) / 2;

      if (progress < 1) {
        // Batch all node updates in a single setNodes call
        setNodes(prevNodes => {
          const updated = prevNodes.map(node => {
            if (startPositions[node.id] && targetPositions[node.id]) {
              const start = startPositions[node.id];
              const target = targetPositions[node.id];
              return {
                ...node,
                x: start.x + (target.x - start.x) * eased,
                y: start.y + (target.y - start.y) * eased
              };
            }
            return node;
          });
          return updated;
        });
        animationFrameId = requestAnimationFrame(animate);
      } else {
        // Final update to ensure exact positions - only once
        setNodes(prevNodes => {
          const updated = prevNodes.map(node => {
            if (startPositions[node.id] && targetPositions[node.id]) {
              const target = targetPositions[node.id];
              return {
                ...node,
                x: target.x,
                y: target.y
              };
            }
            return node;
          });
          return updated;
        });
        setIsAnimating(false);
      }
    };

    animationFrameId = requestAnimationFrame(animate);
    
    // Return cleanup function
    return () => {
      if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
      }
    };
  }, [nodes]);

  /**
   * Handles subgraph highlight without animation.
   * Immediately repositions nodes to show the K5/K3,3 structure.
   * For subdivisions, arranges branch nodes AND their connecting subdivision
   * vertices along straight paths so the true K₅ / K₃,₃ structure is obvious.
   */
  const handleSubgraphHighlight = useCallback((subgraphType) => {
    const solutionInfo = getSolutionInfo();

    if (highlightedSubgraph === subgraphType) {
      // Turning off - restore original positions immediately
      if (originalPositions) {
        setNodes(prevNodes => prevNodes.map(node => {
          if (originalPositions[node.id]) {
            return {
              ...node,
              x: originalPositions[node.id].x,
              y: originalPositions[node.id].y
            };
          }
          return node;
        }));
        setOriginalPositions(null);
      }
      setHighlightedSubgraph(null);
    } else {
      // Turning on - save positions and immediately set to canonical shape
      if (!originalPositions) {
        const saved = {};
        nodes.forEach(node => {
          saved[node.id] = { x: node.x, y: node.y };
        });
        setOriginalPositions(saved);
      }

      setHighlightedSubgraph(subgraphType);

      // Calculate canonical positions based on whether it's a subdivision or exact subgraph
      if (subgraphType === 'k5' && solutionInfo.k5) {
        let canonical = {};

        // 1) Place branch vertices in a clean pentagon
        if (solutionInfo.k5.branchNodes) {
          // It's a subdivision - use branch nodes for canonical positions
          canonical = getCanonicalPositions('k5', solutionInfo.k5.branchNodes) || {};
        } else {
          // It's an exact subgraph
          canonical = getCanonicalPositions('k5', solutionInfo.k5.nodes) || {};
        }

        // 2) For subdivisions, also place the intermediate nodes cleanly along
        // straight segments between the corresponding branch vertices.
        if (solutionInfo.k5.branchNodes && solutionInfo.k5.paths && canonical) {
          const branchPositions = canonical;

          solutionInfo.k5.paths.forEach((path) => {
            if (!Array.isArray(path) || path.length < 2) return;
            const startId = path[0];
            const endId = path[path.length - 1];
            const startPos = branchPositions[startId];
            const endPos = branchPositions[endId];
            if (!startPos || !endPos) return;

            const segments = path.length - 1;
            const dx = endPos.x - startPos.x;
            const dy = endPos.y - startPos.y;

            // Place internal subdivision vertices evenly spaced along the segment.
            for (let i = 1; i < path.length - 1; i++) {
              const nodeId = path[i];
              // t in (0,1) along the segment
              const t = i / segments;
              // Only overwrite if we haven't already positioned this subdivision node
              if (!canonical[nodeId]) {
                canonical[nodeId] = {
                  x: startPos.x + dx * t,
                  y: startPos.y + dy * t
                };
              }
            }
          });
        }

        if (canonical && Object.keys(canonical).length > 0) {
          // Set positions immediately without animation
          setNodes(prevNodes => prevNodes.map(node => {
            if (canonical[node.id]) {
              return {
                ...node,
                x: canonical[node.id].x,
                y: canonical[node.id].y
              };
            }
            return node;
          }));
        }
      } else if (subgraphType === 'k33' && solutionInfo.k33) {
        let canonical = {};

        if (solutionInfo.k33.branchNodes) {
          // It's a subdivision - use branch nodes for canonical positions
          const branchNodes = [...solutionInfo.k33.branchNodes.setA, ...solutionInfo.k33.branchNodes.setB];
          canonical = getCanonicalPositions('k33', branchNodes, {
            setA: solutionInfo.k33.branchNodes.setA,
            setB: solutionInfo.k33.branchNodes.setB
          }) || {};
        } else {
          // It's an exact subgraph
          canonical = getCanonicalPositions('k33', solutionInfo.k33.nodes, solutionInfo.k33) || {};
        }

        // For K₃,₃ subdivisions, do the same: walk each A–B path and line up the
        // internal subdivision vertices between the two branch vertices.
        if (solutionInfo.k33.branchNodes && solutionInfo.k33.paths && canonical) {
          const branchPositions = canonical;

          solutionInfo.k33.paths.forEach((path) => {
            if (!Array.isArray(path) || path.length < 2) return;
            const startId = path[0];
            const endId = path[path.length - 1];
            const startPos = branchPositions[startId];
            const endPos = branchPositions[endId];
            if (!startPos || !endPos) return;

            const segments = path.length - 1;
            const dx = endPos.x - startPos.x;
            const dy = endPos.y - startPos.y;

            for (let i = 1; i < path.length - 1; i++) {
              const nodeId = path[i];
              const t = i / segments;
              if (!canonical[nodeId]) {
                canonical[nodeId] = {
                  x: startPos.x + dx * t,
                  y: startPos.y + dy * t
                };
              }
            }
          });
        }

        if (canonical && Object.keys(canonical).length > 0) {
          // Set positions immediately without animation
          setNodes(prevNodes => prevNodes.map(node => {
            if (canonical[node.id]) {
              return {
                ...node,
                x: canonical[node.id].x,
                y: canonical[node.id].y
              };
            }
            return node;
          }));
        }
      }
    }
  }, [highlightedSubgraph, originalPositions, nodes, getSolutionInfo, getCanonicalPositions]);

  /**
   * Clears highlight and restores positions.
   */
  const clearHighlight = useCallback(() => {
    if (originalPositions) {
      setNodes(prevNodes => prevNodes.map(node => {
        if (originalPositions[node.id]) {
          return {
            ...node,
            x: originalPositions[node.id].x,
            y: originalPositions[node.id].y
          };
        }
        return node;
      }));
      setOriginalPositions(null);
    }
    setHighlightedSubgraph(null);
  }, [originalPositions]);

  // Derived state: is the current embedding planar?
  const isPlanar = metrics.total === 0;

  // --------------------------------------------------------------------------
  // NETWORKX PLANAR LAYOUT - For "Untangle" feature using proper algorithm
  // --------------------------------------------------------------------------

  const [isSolving, setIsSolving] = useState(false);
  const [pyodideStatus, setPyodideStatus] = useState('idle'); // 'idle' | 'loading' | 'ready' | 'error'
  const simulationRef = useRef(null);

  /**
   * Starts the untangle process using NetworkX's planar layout.
   * This uses proper graph algorithms instead of force-directed heuristics.
   */
  const startSolving = useCallback(async () => {
    if (isSolving) return;
    setIsSolving(true);
    setPyodideStatus('loading');

    try {
      // Use NetworkX to compute planar layout
      const result = await computePlanarLayout(nodes, edges);
      setPyodideStatus('ready');

      if (result.isPlanar && result.positions) {
        // Animate to the planar positions
        const targetPositions = nodes.map(node => ({
          id: node.id,
          x: snapToGrid(result.positions[node.id]?.x ?? node.x),
          y: snapToGrid(result.positions[node.id]?.y ?? node.y)
        }));

        // Smooth animation to planar layout
        const startPositions = nodes.map(n => ({ id: n.id, x: n.x, y: n.y }));
        const duration = 1000; // 1 second animation
        const startTime = performance.now();

        const animate = (currentTime) => {
          const elapsed = currentTime - startTime;
          const progress = Math.min(elapsed / duration, 1);
          // Ease out cubic
          const eased = 1 - Math.pow(1 - progress, 3);

          setNodes(startPositions.map((start, i) => {
            const target = targetPositions[i];
            return {
              ...nodes[i],
              x: start.x + (target.x - start.x) * eased,
              y: start.y + (target.y - start.y) * eased
            };
          }));

          if (progress < 1) {
            simulationRef.current = requestAnimationFrame(animate);
          } else {
            // Final snap to grid and analyze
            setNodes(targetPositions.map((pos, i) => ({
              ...nodes[i],
              x: pos.x,
              y: pos.y
            })));
            setIsSolving(false);
            simulationRef.current = null;
            // Re-analyze after animation completes
            setTimeout(() => {
              analyzeGraph(targetPositions.map((pos, i) => ({
                ...nodes[i],
                x: pos.x,
                y: pos.y
              })), edges);
            }, 50);
          }
        };

        simulationRef.current = requestAnimationFrame(animate);
      } else {
        // Graph is not planar - NetworkX confirmed it
        setIsSolving(false);
        setPyodideStatus('ready');
      }
    } catch (error) {
      console.error('Planar layout error:', error);
      setPyodideStatus('error');
      setIsSolving(false);
    }
  }, [isSolving, nodes, edges, analyzeGraph]);

  /**
   * Stops the untangle simulation.
   */
  const stopSolving = useCallback(() => {
    if (simulationRef.current) {
      cancelAnimationFrame(simulationRef.current);
      simulationRef.current = null;
    }
    setIsSolving(false);

    // Snap to grid after solving for cleaner look
    setNodes(prev => {
      const snapped = prev.map(node => ({
        ...node,
        x: snapToGrid(node.x),
        y: snapToGrid(node.y)
      }));
      analyzeGraph(snapped, edges); // Re-analyze after snap
      return snapped;
    });
  }, [analyzeGraph, edges]);

  // Clean up simulation on unmount
  useEffect(() => {
    return () => {
      if (simulationRef.current) {
        cancelAnimationFrame(simulationRef.current);
      }
    };
  }, []);

  // --------------------------------------------------------------------------
  // GLOBAL LOADING OVERLAY - For heavy non-planar graph analysis
  // --------------------------------------------------------------------------

  useEffect(() => {
    // Show the overlay ONLY when:
    // 1) The user has requested the solution,
    // 2) We've already determined the graph is non-planar, and
    // 3) We're currently doing heavy analysis (K₅/K₃,₃ search, fixes, etc.).
    const isConfirmedNonPlanar =
      networkxPlanarity.checked && networkxPlanarity.isPlanar === false;

    const isBusyWithAnalysis =
      networkxPlanarity.checking ||
      checkingMinimumSet ||
      (isSolving && pyodideStatus === 'loading');

    setShowGlobalLoading(showSolution && isConfirmedNonPlanar && isBusyWithAnalysis);
  }, [showSolution, networkxPlanarity, checkingMinimumSet, isSolving, pyodideStatus]);

  // --------------------------------------------------------------------------
  // NODE/EDGE REMOVAL FOR PLANARITY - Find minimum removal sets
  // --------------------------------------------------------------------------

  /**
   * Computes minimum removal sets when graph is non-planar.
   */
  useEffect(() => {
    const computeMinimumRemovalSets = async () => {
      const solutionInfo = getSolutionInfo();
      
      // Only compute for non-planar graphs, and only for reasonably small ones
      if (
        solutionInfo.solvable ||
        nodes.length < 3 ||
        edges.length < 2 ||
        nodes.length > MAX_NODES_FOR_MIN_REMOVAL ||
        edges.length > MAX_EDGES_FOR_MIN_REMOVAL
      ) {
        setMinimumRemovalSet(null);
        setMinimumEdgeRemovalSet(null);
        setRemovalComparison(null);
        return;
      }

      setCheckingMinimumSet(true);
      
      try {
        // Compute both node and edge removal sets in parallel
        const [nodeRemoval, edgeRemoval] = await Promise.all([
          findMinimumNodeRemovalSet(nodes, edges),
          findMinimumEdgeRemovalSet(nodes, edges)
        ]);

        // Generate explanations
        let nodeExplanation = '';
        let edgeExplanation = '';
        
        if (nodeRemoval.nodes.length > 0) {
          nodeExplanation = generateRemovalExplanation(nodeRemoval.nodes, nodes, edges, solutionInfo.k5, solutionInfo.k33);
        }
        
        if (edgeRemoval.edges.length > 0) {
          edgeExplanation = generateEdgeRemovalExplanation(edgeRemoval.edges, nodes, edges, solutionInfo.k5, solutionInfo.k33);
        }

        // Compare strategies
        const comparison = compareRemovalStrategies(nodeRemoval, edgeRemoval);

        setMinimumRemovalSet(nodeRemoval.nodes.length > 0 ? { ...nodeRemoval, explanation: nodeExplanation } : null);
        setMinimumEdgeRemovalSet(edgeRemoval.edges.length > 0 ? { ...edgeRemoval, explanation: edgeExplanation } : null);
        setRemovalComparison(comparison);
      } catch (error) {
        console.error('Error computing minimum removal sets:', error);
        setMinimumRemovalSet(null);
        setMinimumEdgeRemovalSet(null);
        setRemovalComparison(null);
      } finally {
        setCheckingMinimumSet(false);
      }
    };

    // Debounce to avoid excessive computation
    const timer = setTimeout(computeMinimumRemovalSets, 500);
    return () => clearTimeout(timer);
  }, [nodes, edges, getSolutionInfo]);

  /**
   * Removes the minimum set of nodes to make the graph planar and animates to a planar layout.
   */
  const handleRemoveNodeForPlanarity = useCallback(async () => {
    if (isSolving) return;

    // Use minimum removal set
    if (!minimumRemovalSet || minimumRemovalSet.nodes.length === 0) {
      console.warn('No minimum removal set found');
      return;
    }

    const nodesToRemove = minimumRemovalSet.nodes;
    setIsSolving(true);
    setPyodideStatus('loading');

    try {
      // Remove all nodes in the set and all edges connected to them
      const remainingNodes = nodes.filter(n => !nodesToRemove.includes(n.id));
      const remainingEdges = edges.filter(
        e => !nodesToRemove.includes(e.source) && !nodesToRemove.includes(e.target)
      );

      // Check if graph becomes too small
      if (remainingNodes.length < 3) {
        setIsSolving(false);
        setPyodideStatus('ready');
        alert('Removing this node would make the graph too small (less than 3 nodes).');
        return;
      }

      // Renumber nodes to sequential IDs (0, 1, 2, ...)
      const nodeIdMap = new Map();
      const renumberedNodes = remainingNodes.map((node, index) => {
        nodeIdMap.set(node.id, index);
        return { ...node, id: index };
      });

      // Update edge references to match new node IDs
      const renumberedEdges = remainingEdges.map(edge => ({
        source: nodeIdMap.get(edge.source),
        target: nodeIdMap.get(edge.target)
      }));

      // Compute planar layout for the new graph
      const result = await computePlanarLayout(renumberedNodes, renumberedEdges);
      setPyodideStatus('ready');

      // Determine target positions
      let targetPositions;
      if (result.isPlanar && result.positions) {
        // Use computed planar layout
        targetPositions = renumberedNodes.map(node => ({
          id: node.id,
          x: snapToGrid(result.positions[node.id]?.x ?? node.x),
          y: snapToGrid(result.positions[node.id]?.y ?? node.y)
        }));
      } else {
        // Fall back to circular layout
        const circularPositions = generateCircularLayout(renumberedNodes.length);
        targetPositions = renumberedNodes.map((node, index) => ({
          ...node,
          x: snapToGrid(circularPositions[index].x),
          y: snapToGrid(circularPositions[index].y)
        }));
      }

      // Update edges immediately (node removed, edges need to be updated)
      setEdges(renumberedEdges);

      // Update nodes with renumbered IDs but keep their current positions for animation
      const startPositions = renumberedNodes.map(n => ({ id: n.id, x: n.x, y: n.y }));
      setNodes(startPositions);

      // Animate nodes to their new positions
      const duration = 1000; // 1 second animation
      const startTime = performance.now();

      const animate = (currentTime) => {
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);
        // Ease out cubic
        const eased = 1 - Math.pow(1 - progress, 3);

        setNodes(startPositions.map((start, i) => {
          const target = targetPositions[i];
          return {
            ...targetPositions[i],
            x: start.x + (target.x - start.x) * eased,
            y: start.y + (target.y - start.y) * eased
          };
        }));

        if (progress < 1) {
          simulationRef.current = requestAnimationFrame(animate);
        } else {
          // Final snap to grid and analyze
          setNodes(targetPositions);
          setIsSolving(false);
          simulationRef.current = null;
          // Re-analyze after animation completes
          setTimeout(() => {
            analyzeGraph(targetPositions, renumberedEdges);
          }, 50);
        }
      };

      // Start animation on next frame to ensure state updates are applied
      requestAnimationFrame(() => {
        simulationRef.current = requestAnimationFrame(animate);
      });
    } catch (error) {
      console.error('Node removal error:', error);
      setPyodideStatus('error');
      setIsSolving(false);
    }
  }, [nodes, edges, isSolving, minimumRemovalSet, analyzeGraph]);

  /**
   * Removes the minimum set of edges to make the graph planar and animates to a planar layout.
   */
  const handleRemoveEdgesForPlanarity = useCallback(async () => {
    if (isSolving) return;

    // Use minimum edge removal set
    if (!minimumEdgeRemovalSet || minimumEdgeRemovalSet.edges.length === 0) {
      console.warn('No minimum edge removal set found');
      return;
    }

    const edgesToRemove = minimumEdgeRemovalSet.edges;
    const edgeKeys = new Set(edgesToRemove.map(e => 
      `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`
    ));

    setIsSolving(true);
    setPyodideStatus('loading');

    try {
      // Remove the edges
      const remainingEdges = edges.filter(e => {
        const key = `${Math.min(e.source, e.target)}-${Math.max(e.source, e.target)}`;
        return !edgeKeys.has(key);
      });

      // Check if too few edges remain
      if (remainingEdges.length < 2) {
        setIsSolving(false);
        setPyodideStatus('ready');
        alert(`Removing ${edgesToRemove.length} edge${edgesToRemove.length > 1 ? 's' : ''} would leave too few edges.`);
        return;
      }

      // Compute planar layout for the new graph
      const result = await computePlanarLayout(nodes, remainingEdges);
      setPyodideStatus('ready');

      // Determine target positions
      let targetPositions;
      if (result.isPlanar && result.positions) {
        // Use computed planar layout
        targetPositions = nodes.map(node => ({
          id: node.id,
          x: snapToGrid(result.positions[node.id]?.x ?? node.x),
          y: snapToGrid(result.positions[node.id]?.y ?? node.y)
        }));
      } else {
        // Fall back to circular layout
        const circularPositions = generateCircularLayout(nodes.length);
        targetPositions = nodes.map((node, index) => ({
          ...node,
          x: snapToGrid(circularPositions[index].x),
          y: snapToGrid(circularPositions[index].y)
        }));
      }

      // Update edges immediately
      setEdges(remainingEdges);

      // Update nodes with new positions for animation
      const startPositions = nodes.map(n => ({ id: n.id, x: n.x, y: n.y }));
      setNodes(startPositions);

      // Animate nodes to their new positions
      const duration = 1000; // 1 second animation
      const startTime = performance.now();

      const animate = (currentTime) => {
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);
        // Ease out cubic
        const eased = 1 - Math.pow(1 - progress, 3);

        setNodes(startPositions.map((start, i) => {
          const target = targetPositions[i];
          return {
            ...targetPositions[i],
            x: start.x + (target.x - start.x) * eased,
            y: start.y + (target.y - start.y) * eased
          };
        }));

        if (progress < 1) {
          simulationRef.current = requestAnimationFrame(animate);
        } else {
          // Final snap to grid and analyze
          setNodes(targetPositions);
          setIsSolving(false);
          simulationRef.current = null;
          // Re-analyze after animation completes
          setTimeout(() => {
            analyzeGraph(targetPositions, remainingEdges);
          }, 50);
        }
      };

      // Start animation on next frame to ensure state updates are applied
      requestAnimationFrame(() => {
        simulationRef.current = requestAnimationFrame(animate);
      });
    } catch (error) {
      console.error('Edge removal error:', error);
      setPyodideStatus('error');
      setIsSolving(false);
    }
  }, [nodes, edges, isSolving, minimumEdgeRemovalSet, analyzeGraph]);

  // --------------------------------------------------------------------------
  // RENDER - JSX structure of the application
  // --------------------------------------------------------------------------

  return (
    <div className="app">
      {showGlobalLoading && (
        <div className="app-loading-overlay">
          <div className="app-loading-inner">
            <div className="loading-spinner" />
            <div className="app-loading-message">
              Analyzing non-planar structure and searching for obstructions
            </div>
            <div className="app-loading-subtext">
              We are checking this graph\'s structure (including possible K₅ / K₃,₃ witnesses and fixes). This may take a moment, but is carefully limited so the site stays responsive.
            </div>
          </div>
        </div>
      )}
      {/* ================================================================== */}
      {/* HEADER - Logo and title */}
      {/* ================================================================== */}
      <header className="header">
        <div className="header-content">
          <div className="logo-section">
            {/* Triangle graph logo */}
            <div className="logo-icon">
              <svg viewBox="0 0 40 40" width="40" height="40">
                <circle cx="10" cy="10" r="4" fill="currentColor" />
                <circle cx="30" cy="10" r="4" fill="currentColor" />
                <circle cx="20" cy="30" r="4" fill="currentColor" />
                <line x1="10" y1="10" x2="30" y2="10" stroke="currentColor" strokeWidth="2" />
                <line x1="10" y1="10" x2="20" y2="30" stroke="currentColor" strokeWidth="2" />
                <line x1="30" y1="10" x2="20" y2="30" stroke="currentColor" strokeWidth="2" />
              </svg>
            </div>
            <div>
              <h1>Planarity Checker</h1>
            </div>
          </div>
          <button
            className="header-cta-button"
            type="button"
            onClick={() => {
              const base = window.location.origin;
              window.open(`${base}/planarcheckergame (1).html`, "_blank", "noopener,noreferrer");
            }}
          >
            Test your skills &rarr;
          </button>
        </div>
      </header>

      <main className="main-content">
        {/* ================================================================== */}
        {/* CONTROLS BAR - Generation, challenges, and file upload */}
        {/* ================================================================== */}
        <section className="controls-bar">
          <div className="controls-left">
            {/* Quick random generation */}
            <button
              onClick={() => { generateFullyRandom(); setCurrentChallenge(null); setShowSolution(false); setUploadedFileName(null); setHighlightedSubgraph(null); }}
              className="btn btn-primary"
              style={{ display: 'flex', alignItems: 'center', gap: '0.4rem' }}
            >
              <FaRandom /> Random
            </button>

            {/* Custom node count input */}
            <div className="input-with-button">
              <input
                type="number"
                min="2"
                max="1000"
                value={nodeCountInput}
                onChange={(e) => setNodeCountInput(e.target.value)}
                onKeyDown={(e) => { if (e.key === 'Enter') { generateWithNodeCount(); setCurrentChallenge(null); setShowSolution(false); setUploadedFileName(null); setHighlightedSubgraph(null); } }}
                placeholder="Nodes"
              />
              <button onClick={() => { generateWithNodeCount(); setCurrentChallenge(null); setShowSolution(false); setUploadedFileName(null); setHighlightedSubgraph(null); }} className="btn btn-secondary">
                Generate
              </button>
            </div>

            {/* File upload */}
            <input
              type="file"
              ref={fileInputRef}
              onChange={handleFileSelect}
              accept=".json,.graphml,.xml"
              style={{ display: 'none' }}
            />
            <button
              className="btn btn-secondary"
              onClick={handleBrowseClick}
              style={{ display: 'flex', alignItems: 'center', gap: '0.35rem' }}
            >
              <FaUpload /> Upload
            </button>
          </div>

          <div className="controls-right">
            {/* Grid style selector */}
            <select
              value={showGrid ? gridStyle : 'off'}
              onChange={(e) => {
                const value = e.target.value;
                if (value === 'off') {
                  setShowGrid(false);
                } else {
                  setShowGrid(true);
                  setGridStyle(value);
                }
              }}
              className="btn btn-secondary"
              style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', cursor: 'pointer' }}
              title="Grid style (G)"
            >
              <option value="off">No Grid</option>
              <option value="dots">Dots</option>
              <option value="lines">Lines</option>
              <option value="both">Both</option>
            </select>
            {/* Example graphs selector */}
            <select
              value={loadedExampleId || ""}
              onChange={(e) => {
                const graphId = e.target.value;
                if (graphId) {
                  // Load the graph onto the canvas
                  loadExampleGraph(graphId);
                  // Show the analysis/explanation
                  handleAnalyzeGraph(graphId);
                } else {
                  // If "Example Graphs" is selected, clear the loaded graph
                  setLoadedExampleId(null);
                }
              }}
              className="btn btn-secondary"
              style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', cursor: 'pointer', maxWidth: '140px' }}
              title="Load example graph"
            >
              <option value="">Examples</option>
              {importedGraphs.map((graph) => (
                <option key={graph.id} value={graph.id}>{graph.name}</option>
              ))}
            </select>
            {/* Clear graph */}
            {nodes.length > 0 && (
              <button
                className="btn btn-secondary"
                onClick={clearGraph}
                title="Clear graph (Ctrl+Z or Delete)"
                style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '0.35rem' }}
              >
                <FaTrash /> Clear
              </button>
            )}
            {/* Export */}
            {nodes.length > 0 && (
              <div style={{ position: 'relative', display: 'inline-block' }}>
                <button
                  className="btn btn-secondary"
                  onClick={() => setShowExportMenu(!showExportMenu)}
                  onBlur={() => setTimeout(() => setShowExportMenu(false), 200)}
                  title="Export graph"
                  style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '0.35rem' }}
                >
                  <FaDownload /> Export
                </button>
                {showExportMenu && (
                  <div className="export-dropdown">
                    <button 
                      onClick={() => { exportGraph(); setShowExportMenu(false); }} 
                      className="export-option"
                    >
                      📄 JSON
                    </button>
                    <button 
                      onClick={() => { exportAsImage(); setShowExportMenu(false); }} 
                      className="export-option"
                    >
                      🖼️ PNG
                    </button>
                  </div>
                )}
              </div>
            )}
            {/* 3D View toggle */}
            <button
              className={`btn btn-3d ${show3DView ? 'active' : ''}`}
              onClick={() => setShow3DView(!show3DView)}
              disabled={nodes.length === 0}
              style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '0.35rem' }}
            >
              <FaGlobe /> {show3DView ? 'Hide 3D' : '3D'}
            </button>
            {/* Solution toggle */}
            <button
              className="btn btn-solution"
              onClick={() => {
                setShowSolution(!showSolution);
                if (showSolution) clearHighlight(); // Clear highlight when hiding
              }}
              style={{ padding: '0.5rem 0.75rem', fontSize: '0.8rem', display: 'flex', alignItems: 'center', gap: '0.35rem' }}
            >
              <FaLightbulb /> {showSolution ? 'Hide' : 'Solution'}
            </button>
          </div>
        </section>

        <div className="top-layout-row">
        {/* Challenge banner */}
        {currentChallenge && (
          <div className="challenge-banner-top">
            <span className="challenge-label">Challenge: {EXAMPLE_GRAPHS[currentChallenge].name}</span>
            <span className="challenge-status">
              {isPlanar ? '🎉 Impossible! You found a planar embedding!' : `${metrics.total} crossings remaining`}
            </span>
          </div>
        )}

        {/* Solution content */}
        {showSolution && (() => {
          const solutionInfo = getSolutionInfo();
          return (
            <div className={`solution-content-top ${solutionInfo.solvable ? 'solvable' : 'impossible'}`}>
              <div className="solution-verdict">
                <span className="verdict-icon">
                  {solutionInfo.loading ? '⏳' : solutionInfo.solvable ? '✓' : '✗'}
                </span>
                <span className="verdict-text">
                  {solutionInfo.loading ? 'Checking...' : solutionInfo.solvable ? 'Planar (Solvable)' : 'Non-Planar (Impossible)'}
                </span>
                {planarityCheckTime !== null && !solutionInfo.loading && (
                  <div className="timer-box">
                    <span className="timer-label">Time:</span>
                    <span className="timer-value">
                      {planarityCheckTime < 1 
                        ? `${(planarityCheckTime * 1000).toFixed(1)} μs`
                        : planarityCheckTime < 1000
                          ? `${planarityCheckTime.toFixed(2)} ms`
                          : `${(planarityCheckTime / 1000).toFixed(3)} s`}
                    </span>
                  </div>
                )}
              </div>
              <p className="solution-reason">{solutionInfo.reason}</p>
              
              {/* Thought Process Button */}
              {solutionInfo.thoughtProcess && solutionInfo.thoughtProcess.length > 0 && (
                <div style={{ marginTop: '1rem' }}>
                  <button
                    className="btn btn-thought-process"
                    onClick={() => setShowThoughtProcess(!showThoughtProcess)}
                  >
                    {showThoughtProcess ? (
                      <>
                        <FaChevronUp style={{ marginRight: '0.4rem' }} />
                        <span>Hide Thought Process</span>
                      </>
                    ) : (
                      <>
                        <FaBrain style={{ marginRight: '0.4rem' }} />
                        <span>Show Thought Process</span>
                      </>
                    )}
                  </button>
                </div>
              )}
              
              {/* Thought Process Display */}
              {showThoughtProcess && solutionInfo.thoughtProcess && solutionInfo.thoughtProcess.length > 0 && (
                <div className="thought-process-container">
                  <h4 style={{ marginBottom: '0.75rem', fontSize: '0.9rem', fontWeight: 600, color: 'var(--accent-primary)' }}>
                    Algorithm Thought Process:
                  </h4>
                  {solutionInfo.thoughtProcess.map((step, index) => (
                    <div key={index} className="thought-process-step">
                      <span className="thought-process-step-number">{step.step}.</span>
                      <div className="thought-process-step-content">
                        <div className="thought-process-step-description">{step.description}</div>
                        <div className="thought-process-step-result">{step.result}</div>
                      </div>
                    </div>
                  ))}
                </div>
              )}

              {/* Untangle Button - Only for solvable graphs */}
              {solutionInfo.solvable && (
                <div className="untangle-section" style={{ marginTop: '1rem', marginBottom: '1rem' }}>
                  <button
                    className={`btn btn-primary ${isSolving ? 'active' : ''}`}
                    onClick={isSolving ? stopSolving : startSolving}
                    disabled={isSolving}
                  >
                    {isSolving && pyodideStatus === 'loading'
                      ? (
                        <>
                          <span style={{ marginRight: '0.4rem' }}>⏳</span>
                          Loading NetworkX...
                        </>
                      )
                      : isSolving
                        ? (
                          <>
                            <FaMagic style={{ marginRight: '0.4rem' }} />
                            Computing planar layout...
                          </>
                        )
                        : (
                          <>
                            <FaMagic style={{ marginRight: '0.4rem' }} />
                            Untangle Graph
                          </>
                        )}
                  </button>
                </div>
              )}

              {/* Famous Graph Info */}
              {solutionInfo.famousGraph && (
                <div className="famous-graph-info">
                  <strong>Famous Graph Identified:</strong> {solutionInfo.famousGraph.name}
                </div>
              )}

              {/* Fix Suggestions - Only show for solvable (planar) graphs with crossings */}
              {solutionInfo.solvable && solutionInfo.fixes && (
                <div className="fix-suggestions">
                  <p className="fix-hint">
                    Found {solutionInfo.fixes.nodes.size} nodes and {solutionInfo.fixes.edges.size} edges that can be removed to fix planarity issues.
                  </p>
                  <button
                    className={`btn btn-fix ${showFixes ? 'active' : ''}`}
                    onClick={() => {
                      setShowFixes(!showFixes);
                      if (!showFixes) setHighlightedSubgraph(null); // Clear other highlights
                    }}
                  >
                    {showFixes ? 'Hide Fixes' : '✨ Show Fix Suggestions'}
                  </button>
                </div>
              )}

              {/* Subgraph highlight buttons */}
              {(solutionInfo.k5 || solutionInfo.k33) && (
                <div className="subgraph-buttons">
                  <span className="subgraph-label">Contains:</span>
                  {solutionInfo.k5 && (
                    <button
                      className={`btn btn-subgraph ${highlightedSubgraph === 'k5' ? 'active' : ''}`}
                      onClick={() => {
                        handleSubgraphHighlight('k5');
                        setShowFixes(false);
                      }}
                    >
                      K₅ {solutionInfo.k5.branchNodes ? '(subdivision)' : '(subgraph)'} ({solutionInfo.k5.nodes.length} nodes)
                    </button>
                  )}
                  {solutionInfo.k33 && (
                    <button
                      className={`btn btn-subgraph ${highlightedSubgraph === 'k33' ? 'active' : ''}`}
                      onClick={() => {
                        handleSubgraphHighlight('k33');
                        setShowFixes(false);
                      }}
                    >
                      K₃,₃ {solutionInfo.k33.branchNodes ? '(subdivision)' : '(subgraph)'} ({solutionInfo.k33.nodes.length} nodes)
                    </button>
                  )}
                  {(highlightedSubgraph || showFixes) && (
                    <button
                      className="btn btn-clear-highlight"
                      onClick={() => {
                        clearHighlight();
                        setShowFixes(false);
                      }}
                    >
                      Clear
                    </button>
                  )}
                </div>
              )}

              {/* Removal Strategy Comparison - Show both node and edge removal options */}
              {!solutionInfo.solvable && (minimumRemovalSet || minimumEdgeRemovalSet) && (
                <div className="untangle-section" style={{ marginTop: '1rem', marginBottom: '1rem' }}>
                  {/* Comparison header */}
                  {removalComparison && (
                    <div style={{ 
                      padding: '0.75rem', 
                      background: removalComparison.better === 'nodes' 
                        ? 'rgba(239, 68, 68, 0.1)' 
                        : removalComparison.better === 'edges'
                        ? 'rgba(16, 185, 129, 0.1)'
                        : 'rgba(100, 116, 139, 0.1)',
                      border: `1px solid ${removalComparison.better === 'nodes' 
                        ? 'rgba(239, 68, 68, 0.3)' 
                        : removalComparison.better === 'edges'
                        ? 'rgba(16, 185, 129, 0.3)'
                        : 'rgba(100, 116, 139, 0.3)'}`, 
                      borderRadius: '8px',
                      marginBottom: '1rem'
                    }}>
                      <p style={{ margin: '0 0 0.5rem 0', fontWeight: 600, 
                        color: removalComparison.better === 'nodes' 
                          ? 'var(--danger)' 
                          : removalComparison.better === 'edges'
                          ? 'var(--success)'
                          : 'var(--text-secondary)' }}>
                        {removalComparison.better === 'nodes' && '✓ Recommended: Remove Nodes'}
                        {removalComparison.better === 'edges' && '✓ Recommended: Remove Edges'}
                        {removalComparison.better === 'equal' && 'Both Options Available'}
                        {removalComparison.better === 'none' && 'No Viable Solution'}
                      </p>
                      <p style={{ margin: 0, fontSize: '0.85rem', color: 'var(--text-secondary)', lineHeight: '1.5' }}>
                        {removalComparison.reason}
                      </p>
                    </div>
                  )}

                  {/* Node Removal Option */}
                  {minimumRemovalSet && minimumRemovalSet.nodes.length > 0 && (
                    <div style={{ 
                      padding: '0.75rem', 
                      background: removalComparison?.better === 'nodes' 
                        ? 'rgba(239, 68, 68, 0.05)' 
                        : 'rgba(59, 130, 246, 0.1)', 
                      border: `1px solid ${removalComparison?.better === 'nodes' 
                        ? 'rgba(239, 68, 68, 0.3)' 
                        : 'rgba(59, 130, 246, 0.3)'}`, 
                      borderRadius: '8px',
                      marginBottom: '1rem'
                    }}>
                      <p style={{ margin: '0 0 0.5rem 0', fontWeight: 600, color: 'var(--accent-primary)' }}>
                        Option 1: Remove {minimumRemovalSet.count} Node{minimumRemovalSet.count > 1 ? 's' : ''}
                      </p>
                      <p style={{ margin: '0 0 0.5rem 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                        <strong>Nodes to remove:</strong> {minimumRemovalSet.nodes.join(', ')}
                      </p>
                      <p style={{ margin: '0 0 0.75rem 0', fontSize: '0.85rem', color: 'var(--text-secondary)', lineHeight: '1.5' }}>
                        {minimumRemovalSet.explanation}
                      </p>
                      <button
                        className={`btn ${removalComparison?.better === 'nodes' ? 'btn-primary' : 'btn-secondary'}`}
                        onClick={handleRemoveNodeForPlanarity}
                        disabled={isSolving || checkingMinimumSet}
                        style={{ width: '100%' }}
                      >
                        {checkingMinimumSet ? (
                          <>
                            <FaBrain style={{ marginRight: '0.4rem' }} />
                            <span>Finding minimum removal set...</span>
                          </>
                        ) : isSolving && pyodideStatus === 'loading' ? (
                          '⏳ Computing...'
                        ) : isSolving ? (
                          `✨ Removing ${minimumRemovalSet.count} node${minimumRemovalSet.count > 1 ? 's' : ''}...`
                        ) : (
                          <>
                            <FaTrash style={{ marginRight: '0.4rem' }} />
                            <span>
                              Remove {minimumRemovalSet.count} Node{minimumRemovalSet.count > 1 ? 's' : ''}
                            </span>
                          </>
                        )}
                      </button>
                    </div>
                  )}

                  {/* Edge Removal Option */}
                  {minimumEdgeRemovalSet && minimumEdgeRemovalSet.edges.length > 0 && (
                    <div style={{ 
                      padding: '0.75rem', 
                      background: removalComparison?.better === 'edges' 
                        ? 'rgba(16, 185, 129, 0.05)' 
                        : 'rgba(59, 130, 246, 0.1)', 
                      border: `1px solid ${removalComparison?.better === 'edges' 
                        ? 'rgba(16, 185, 129, 0.3)' 
                        : 'rgba(59, 130, 246, 0.3)'}`, 
                      borderRadius: '8px',
                      marginBottom: '1rem'
                    }}>
                      <p style={{ margin: '0 0 0.5rem 0', fontWeight: 600, color: 'var(--accent-primary)' }}>
                        Option 2: Remove {minimumEdgeRemovalSet.count} Edge{minimumEdgeRemovalSet.count > 1 ? 's' : ''}
                      </p>
                      <p style={{ margin: '0 0 0.5rem 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                        <strong>Edges to remove:</strong> {minimumEdgeRemovalSet.edges.map(e => `${e.source}-${e.target}`).join(', ')}
                      </p>
                      <p style={{ margin: '0 0 0.75rem 0', fontSize: '0.85rem', color: 'var(--text-secondary)', lineHeight: '1.5' }}>
                        {minimumEdgeRemovalSet.explanation}
                      </p>
                      <button
                        className={`btn ${removalComparison?.better === 'edges' ? 'btn-primary' : 'btn-secondary'}`}
                        onClick={handleRemoveEdgesForPlanarity}
                        disabled={isSolving || checkingMinimumSet}
                        style={{ width: '100%' }}
                      >
                        {checkingMinimumSet ? (
                          <>
                            <FaBrain style={{ marginRight: '0.4rem' }} />
                            <span>Finding minimum removal set...</span>
                          </>
                        ) : isSolving && pyodideStatus === 'loading' ? (
                          '⏳ Computing...'
                        ) : isSolving ? (
                          `✨ Removing ${minimumEdgeRemovalSet.count} edge${minimumEdgeRemovalSet.count > 1 ? 's' : ''}...`
                        ) : (
                          <>
                            <FaCut style={{ marginRight: '0.4rem' }} />
                            <span>
                              Remove {minimumEdgeRemovalSet.count} Edge{minimumEdgeRemovalSet.count > 1 ? 's' : ''}
                            </span>
                          </>
                        )}
                      </button>
                    </div>
                  )}

                  {checkingMinimumSet && (
                    <p className="fix-hint" style={{ marginTop: '0.5rem', textAlign: 'center' }}>
                      Comparing node removal vs edge removal strategies...
                    </p>
                  )}
                </div>
              )}
            </div>
          );
        })()}

        {/* ================================================================== */}
        {/* 3D VIEW SECTION - Sphere visualization of the graph */}
        {/* ================================================================== */}
        {!loadedExampleId && show3DView && nodes.length > 0 && (
          <section className="view-3d-section">
            <div className="view-3d-header">
              <h3>🌐 3D Sphere View</h3>
              <div className="view-3d-controls">
                <button
                  className="btn btn-small"
                  onClick={() => setSphereRotation({ x: sphereRotation.x - 15, y: sphereRotation.y })}
                >
                  ↑
                </button>
                <button
                  className="btn btn-small"
                  onClick={() => setSphereRotation({ x: sphereRotation.x + 15, y: sphereRotation.y })}
                >
                  ↓
                </button>
                <button
                  className="btn btn-small"
                  onClick={() => setSphereRotation({ x: sphereRotation.x, y: sphereRotation.y - 15 })}
                >
                  ←
                </button>
                <button
                  className="btn btn-small"
                  onClick={() => setSphereRotation({ x: sphereRotation.x, y: sphereRotation.y + 15 })}
                >
                  →
                </button>
                <button
                  className="btn btn-small btn-secondary"
                  onClick={() => setSphereRotation({ x: 0, y: 0 })}
                >
                  Reset
                </button>
                <button
                  className={`btn btn-small ${autoRotate3D ? 'btn-active' : 'btn-secondary'}`}
                  onClick={() => setAutoRotate3D(!autoRotate3D)}
                  title={autoRotate3D ? 'Disable auto-rotate' : 'Enable auto-rotate'}
                >
                  {autoRotate3D ? '⏸ Auto' : '▶ Auto'}
                </button>
              </div>
            </div>
            <div className="view-3d-container">
              <svg
                className={`view-3d-svg ${isDragging3D ? 'dragging' : ''}`}
                viewBox="0 0 400 400"
                onMouseDown={(e) => {
                  setIsDragging3D(true);
                  setLastMousePos({ x: e.clientX, y: e.clientY });
                }}
                onMouseMove={(e) => {
                  if (!isDragging3D) return;
                  const dx = e.clientX - lastMousePos.x;
                  const dy = e.clientY - lastMousePos.y;
                  setSphereRotation(prev => ({
                    x: prev.x + dy * 0.5,
                    y: prev.y + dx * 0.5
                  }));
                  setLastMousePos({ x: e.clientX, y: e.clientY });
                }}
                onMouseUp={() => setIsDragging3D(false)}
                onMouseLeave={() => setIsDragging3D(false)}
              >
                {/* Background sphere wireframe */}
                <circle cx="200" cy="200" r="120" className="sphere-outline" />
                <ellipse cx="200" cy="200" rx="120" ry="40" className="sphere-ring" />
                <ellipse cx="200" cy="200" rx="40" ry="120" className="sphere-ring" />

                {/* Draw edges (back to front) */}
                {get3DEdges().map((edge, idx) => (
                  <line
                    key={`3d-edge-${idx}`}
                    x1={edge.x1}
                    y1={edge.y1}
                    x2={edge.x2}
                    y2={edge.y2}
                    className="edge-3d"
                    style={{ opacity: edge.opacity }}
                  />
                ))}

                {/* Draw nodes (back to front by z-index) */}
                {get3DNodePositions()
                  .sort((a, b) => a.z - b.z)
                  .map((node) => (
                    <g key={`3d-node-${node.id}`}>
                      <circle
                        cx={node.x}
                        cy={node.y}
                        r={4 + node.scale * 3}
                        className="node-3d"
                        style={{ opacity: node.opacity }}
                      />
                      {nodes.length <= 30 && (
                        <text
                          x={node.x}
                          y={node.y + 1}
                          className="node-label-3d"
                          style={{
                            opacity: node.opacity,
                            fontSize: `${6 + node.scale * 2}px`
                          }}
                        >
                          {node.id}
                        </text>
                      )}
                    </g>
                  ))}
              </svg>
              <p className="view-3d-hint">Drag to rotate • {autoRotate3D ? 'Auto-rotating' : 'Auto-rotate off'}</p>
            </div>
          </section>
        )}

        {/* ================================================================== */}
        {/* GRAPH SECTION - Main interactive graph area */}
        {/* ================================================================== */}
        {!loadedExampleId && (
        <section className={`graph-section ${(showSolution || show3DView) ? 'graph-section-right' : ''}`}>
          {/* Graph container */}
          <div
            className={`graph-container ${isDragOver ? 'drag-over' : ''} ${nodes.length > 0 ? 'has-graph' : ''}`}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
          >
            <svg
              className={`graph-svg ${isPanning2D ? 'panning' : ''}`}
              viewBox={getViewBox()}
              ref={svgRef}
              onMouseDown={handleCanvasMouseDown}
              onMouseMove={handleMouseMove}
              onMouseUp={handleMouseUp}
              onMouseLeave={handleMouseUp}
              onTouchStart={handleTouchStart}
              onTouchMove={handleTouchMove}
              onTouchEnd={handleTouchEnd}
              onWheel={handleWheel}
            >
              {/* Grid background - improved with toggle and styles */}
              {showGrid && (
                <>
                  <defs>
                    {(gridStyle === 'lines' || gridStyle === 'both') && (
                      <pattern id="gridLines" width={GRID_SIZE} height={GRID_SIZE} patternUnits="userSpaceOnUse">
                        <path d={`M ${GRID_SIZE} 0 L 0 0 0 ${GRID_SIZE}`} fill="none" stroke="#3b82f6" strokeWidth="0.5" opacity="0.15" />
                      </pattern>
                    )}
                    {/* Major grid lines every 5 cells */}
                    {(gridStyle === 'lines' || gridStyle === 'both') && (
                      <pattern id="gridLinesMajor" width={GRID_SIZE * 5} height={GRID_SIZE * 5} patternUnits="userSpaceOnUse">
                        <path d={`M ${GRID_SIZE * 5} 0 L 0 0 0 ${GRID_SIZE * 5}`} fill="none" stroke="#3b82f6" strokeWidth="0.8" opacity="0.25" />
                      </pattern>
                    )}
                  </defs>
                  {(gridStyle === 'lines' || gridStyle === 'both') && (
                    <>
                      <rect className="grid-background" width="400" height="400" fill="url(#gridLinesMajor)" style={{ cursor: 'grab' }} />
                      <rect className="grid-background" width="400" height="400" fill="url(#gridLines)" style={{ cursor: 'grab' }} />
                    </>
                  )}
                  
                  {/* Grid dots - snap points where nodes can be placed */}
                  {(gridStyle === 'dots' || gridStyle === 'both') && nodes.length < 200 && getGridPoints().map((point, index) => {
                    const isOccupied = nodes.some(n => n.x === point.x && n.y === point.y);
                    return (
                      <circle
                        key={`grid-${index}`}
                        cx={point.x}
                        cy={point.y}
                        r={isOccupied ? 1 : 2}
                        fill={isOccupied ? 'transparent' : '#3b82f6'}
                        opacity={isOccupied ? 0 : 0.25}
                        className="grid-dot"
                      />
                    );
                  })}
                </>
              )}
              {!showGrid && (
                <rect className="grid-background" width="400" height="400" fill="transparent" style={{ cursor: 'grab' }} />
              )}

              {/* Draw edges */}
              {(() => {
                const highlightedEdgeSet = getHighlightedEdges();

                return edges.map((edge, index) => {
                  const startNode = nodes.find(n => n.id === edge.source);
                  const endNode = nodes.find(n => n.id === edge.target);
                  if (!startNode || !endNode) return null;
                  const crossings = metrics.edgeScores[index] || 0;
                  const isBad = crossings > 0;

                  // Check if this edge is part of highlighted subgraph
                  const edgeKey = `${Math.min(edge.source, edge.target)}-${Math.max(edge.source, edge.target)}`;
                  const isHighlighted = highlightedEdgeSet.has(edgeKey);
                  const hasHighlight = highlightedEdgeSet.size > 0;

                  let className = 'edge';
                  if (isHighlighted) className = showFixes ? 'edge highlighted-fix' : 'edge highlighted-subgraph';
                  else if (hasHighlight) className = 'edge dimmed'; // Dim non-highlighted edges
                  else if (isBad) className = 'edge crossing';

                  const isHovered = hoveredEdge === index;
                  return (
                    <line
                      key={`edge-${index}`}
                      x1={startNode.x} y1={startNode.y}
                      x2={endNode.x} y2={endNode.y}
                      className={`${className} ${isHovered ? 'edge-hovered' : ''}`}
                      onMouseEnter={() => setHoveredEdge(index)}
                      onMouseLeave={() => setHoveredEdge(null)}
                      style={{ cursor: 'pointer' }}
                    />
                  );
                });
              })()}

              {/* Draw draggable nodes - size adapts to node count */}
              {(() => {
                // Dynamic node sizing based on count
                const nodeRadius = nodes.length > 200 ? 2 : nodes.length > 50 ? 3 : nodes.length > 20 ? 4 : 5;
                const showLabels = nodes.length <= 50;
                const highlightedNodeSet = getHighlightedNodes(); // Only branch nodes (5 or 6)
                const subdivisionNodeSet = getSubdivisionNodes(); // Intermediate nodes on paths
                const branchNodeSet = getBranchNodes();
                const hasHighlight = highlightedNodeSet.size > 0;

                return nodes.map((node) => {
                  const isBranchNode = branchNodeSet.has(node.id);
                  const isSubdivisionNode = subdivisionNodeSet.has(node.id);
                  const isRelevant = isBranchNode || isSubdivisionNode;
                  const isCurrentlyMoving = currentlyMovingNode === node.id;

                  let nodeClass = 'node';
                  if (isCurrentlyMoving) nodeClass = 'node currently-moving';
                  else if (isBranchNode) nodeClass = 'node highlighted-branch-node';
                  else if (isSubdivisionNode) nodeClass = 'node highlighted-subdivision-node';
                  else if (hasHighlight) nodeClass = 'node dimmed-irrelevant';

                  // Make moving node larger and more visible
                  // Branch nodes stay smaller and fixed when highlighted
                  const displayRadius = isCurrentlyMoving
                    ? nodeRadius + 4
                    : isBranchNode
                      ? nodeRadius  // Keep branch nodes same size as regular nodes
                      : isSubdivisionNode
                        ? nodeRadius + 1
                        : nodeRadius;

                  const isSelected = selectedNode === node.id;
                  const isDragging = draggingNodeRef.current === node.id && isDraggingNode;
                  return (
                    <g
                      key={`node-${node.id}`}
                      onMouseDown={isBranchNode ? undefined : (e) => {
                        handleMouseDown(e, node.id);
                        setSelectedNode(node.id);
                      }}
                      className={isBranchNode ? "" : `draggable-node ${isDragging ? 'dragging' : ''}`}
                      style={isBranchNode ? { cursor: 'default', pointerEvents: 'none' } : { cursor: isDragging ? 'grabbing' : 'grab' }}
                    >
                      <circle
                        cx={node.x}
                        cy={node.y}
                        r={displayRadius}
                        className={nodeClass}
                      />
                      {showLabels && (
                        <text
                          x={node.x}
                          y={node.y + nodeRadius * 0.4}
                          className={`node-label ${hasHighlight && !isRelevant && !isCurrentlyMoving ? 'dimmed-irrelevant' : ''}`}
                          style={{ fontSize: `${Math.max(5, nodeRadius * 0.8)}px` }}
                        >
                          {node.id}
                        </text>
                      )}
                    </g>
                  );
                });
              })()}
            </svg>

            {/* Stats overlay when graph is loaded */}
            {nodes.length > 0 && (
              <div className="graph-stats-overlay">
                <div className="graph-stat">
                  <span className="graph-stat-value">{nodes.length}</span>
                  <span className="graph-stat-label">Nodes</span>
                </div>
                <div className="graph-stat">
                  <span className="graph-stat-value">{edges.length}</span>
                  <span className="graph-stat-label">Edges</span>
                </div>
                <div className="graph-stat">
                  <span className="graph-stat-value">{metrics.total}</span>
                  <span className="graph-stat-label">Crossings</span>
                </div>
                <div className={`graph-stat verdict ${isPlanar ? 'planar' : 'non-planar'}`}>
                  <span className="graph-stat-icon">{isPlanar ? '✓' : '✗'}</span>
                  <span className="graph-stat-label">{isPlanar ? 'Planar' : 'Non-Planar'}</span>
                </div>
              </div>
            )}

            {/* Zoom controls */}
            <div className="zoom-controls">
              <button className="zoom-btn" onClick={zoomOut} title="Zoom Out">−</button>
              <button className="zoom-btn zoom-reset" onClick={resetZoom} title="Reset Zoom">{Math.round(zoom * 100)}%</button>
              <button className="zoom-btn" onClick={zoomIn} title="Zoom In">+</button>
            </div>

            {/* Drop overlay hint - only show when no graph */}
            {nodes.length === 0 && (
              <div className="graph-drop-hint">
                <svg viewBox="0 0 32 32" width="48" height="48" fill="none" stroke="currentColor" strokeWidth="1.5" style={{ marginBottom: '0.5rem' }}>
                  <path d="M16 22V10M16 10L11 15M16 10L21 15" />
                  <path d="M6 22V24C6 25.1 6.9 26 8 26H24C25.1 26 26 25.1 26 24V22" />
                </svg>
                <span style={{ fontSize: '1.1rem', fontWeight: '600', marginBottom: '0.5rem' }}>Get Started</span>
                <span style={{ marginBottom: '1rem' }}>Upload a graph file or generate a new graph</span>
                <div style={{ display: 'flex', gap: '0.75rem', marginTop: '0.5rem' }}>
                  <button
                    className="btn btn-primary"
                    onClick={handleBrowseClick}
                    style={{ padding: '0.5rem 1rem', fontSize: '0.875rem', display: 'flex', alignItems: 'center', gap: '0.4rem' }}
                  >
                    <FaUpload /> Upload File
                  </button>
                  <button
                    className="btn btn-secondary"
                    onClick={() => {
                      generateFullyRandom();
                      setCurrentChallenge(null);
                      setShowSolution(false);
                      setUploadedFileName(null);
                      setHighlightedSubgraph(null);
                    }}
                    style={{ padding: '0.5rem 1rem', fontSize: '0.875rem' }}
                  >
                    🎲 Generate Graph
                  </button>
                </div>
                <span className="drop-hint-sub" style={{ marginTop: '0.75rem' }}>Or drag and drop a file here</span>
              </div>
            )}

            {/* Upload status feedback */}
            {uploadStatus && (
              <div className={`upload-status-overlay ${uploadStatus.type}`}>
                <span className="status-icon">{uploadStatus.type === 'success' ? '✓' : '✗'}</span>
                <span className="status-message">{uploadStatus.message}</span>
              </div>
            )}
          </div>
        </section>
        )}
        </div>

        {/* ================================================================== */}
        {/* ANALYSIS PANEL - Shows when an example graph is analyzed */}
        {/* ================================================================== */}
        {activeAnalysis && (
          <section className="analysis-section">
            <div className="analysis-header">
              <div>
                <h2>{activeAnalysis.name}</h2>
                <p>{activeAnalysis.description}</p>
              </div>
            </div>

            <div className="analysis-content">
              {/* Graph visualization */}
              <div className="analysis-graph">
                <svg viewBox="0 0 400 400" width="100%" height="100%">
                  {/* Draw edges first (below nodes) */}
                  {activeAnalysis.edges.map((edge, index) => {
                    const startNode = activeAnalysis.nodes.find(n => n.id === edge.source);
                    const endNode = activeAnalysis.nodes.find(n => n.id === edge.target);
                    if (!startNode || !endNode) return null;
                    const isCrossing = activeAnalysis.edgeScores[index] > 0;
                    return (
                      <line
                        key={index}
                        x1={startNode.x} y1={startNode.y}
                        x2={endNode.x} y2={endNode.y}
                        className={isCrossing ? 'edge crossing' : 'edge'}
                      />
                    );
                  })}
                  {/* Draw nodes on top */}
                  {activeAnalysis.nodes.map((node) => (
                    <g key={node.id}>
                      <circle cx={node.x} cy={node.y} r="12" className="node" />
                      <text x={node.x} y={node.y + 4} className="node-label">{node.id}</text>
                    </g>
                  ))}
                </svg>
              </div>

              {/* Analysis information panel */}
              <div className="analysis-info">
                {/* Planar/Non-planar verdict */}
                <div className={`analysis-result ${activeAnalysis.isPlanar ? 'planar' : 'non-planar'}`}>
                  <div className="result-icon">{activeAnalysis.isPlanar ? '✓' : '✗'}</div>
                  <div className="result-text">
                    <span className="result-label">{activeAnalysis.isPlanar ? 'Planar' : 'Non-Planar'}</span>
                    <span className="result-crossings">{activeAnalysis.crossings} edge crossings detected</span>
                  </div>
                </div>

                {/* Explanation of why it's planar/non-planar */}
                <div className="analysis-explanation">
                  <h4>Why?</h4>
                  <p>{activeAnalysis.explanation}</p>
                </div>

                {/* Graph statistics */}
                <div className="analysis-stats">
                  <div className="analysis-stat">
                    <span className="stat-value">{activeAnalysis.nodes.length}</span>
                    <span className="stat-label">Vertices</span>
                  </div>
                  <div className="analysis-stat">
                    <span className="stat-value">{activeAnalysis.edges.length}</span>
                    <span className="stat-label">Edges</span>
                  </div>
                  <div className="analysis-stat">
                    <span className="stat-value">{activeAnalysis.crossings}</span>
                    <span className="stat-label">Crossings</span>
                  </div>
                </div>
              </div>
            </div>
          </section>
        )}

      </main>
    </div>
  );
}

export default App;
