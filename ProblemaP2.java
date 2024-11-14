import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;

class ProblemaP2 {

    private static Set<Integer> maxFlowVertices = new HashSet<>();
    private static Map<Integer, Integer> indexToCellId = new HashMap<>();

    public static void main(String[] args) throws Exception {
        long startTime = System.currentTimeMillis();

        ProblemaP2 instancia = new ProblemaP2();
        try (
                InputStreamReader is = new InputStreamReader(System.in);
                BufferedReader br = new BufferedReader(is)) {
            String line = br.readLine();
            int casos = Integer.parseInt(line);

            for (int i = 0; i < casos; i++) {
                line = br.readLine();
                String[] primeraLinea = line.split(" ");
                int n = Integer.parseInt(primeraLinea[0]);
                int d = Integer.parseInt(primeraLinea[1]);
                List<Celula> celulas = new ArrayList<>();
                for (int j = 0; j < n; j++) {
                    line = br.readLine();
                    String[] dataStr = line.split(" ");
                    int id = Integer.parseInt(dataStr[0]);
                    int x = Integer.parseInt(dataStr[1]);
                    int y = Integer.parseInt(dataStr[2]);
                    int tipo = Integer.parseInt(dataStr[3]);
                    Celula celula = instancia.new Celula(id, x, y, tipo);
                    for (int k = 4; k < dataStr.length; k++) {
                        celula.peptidos.add(dataStr[k]);
                    }
                    celulas.add(celula);
                }
                Resultado resultado = resolverCaso(celulas, d);
                System.out.println(resultado.celulaId + " " + resultado.flujoOriginal + " " + resultado.flujoReducido);
            }
        }
    }

    private static class Resultado {
        int celulaId;
        int flujoOriginal;
        int flujoReducido;

        public Resultado(int celulaId, int flujoOriginal, int flujoReducido) {
            this.celulaId = celulaId;
            this.flujoOriginal = flujoOriginal;
            this.flujoReducido = flujoReducido;
        }
    }

    public class Celula {
        int id;
        int posX;
        int posY;
        int tipo;
        List<String> peptidos;

        public Celula(int id, int posX, int posY, int tipo) {
            this.id = id;
            this.posX = posX;
            this.posY = posY;
            this.tipo = tipo;
            this.peptidos = new ArrayList<>();
        }
    }

    private static Resultado resolverCaso(List<Celula> celulas, double d) {
        int n = celulas.size();
        int source = n;
        int sink = n + 1;
        Map<Integer, Integer> cellIdToIndex = new HashMap<>();
        indexToCellId.clear();
        List<List<Edge>> graph = buildGraph(celulas, d, cellIdToIndex, source, sink);

        maxFlowVertices.clear();
        int flujoOriginal = maxFlow(graph, source, sink, -1);

        int mejorCelulaId = -1;
        int menorFlujo = Integer.MAX_VALUE;

        for (int i = 0; i < celulas.size(); i++) {
            Celula celula = celulas.get(i);
            if (celula.tipo == 2 && maxFlowVertices.contains(celula.id)) { // Calculadora
                int bloqueadaIndex = cellIdToIndex.get(celula.id);
                int flujoActual = maxFlow(graph, source, sink, bloqueadaIndex);
                if (flujoActual <= menorFlujo) {
                    menorFlujo = flujoActual;
                    mejorCelulaId = celula.id;
                }
            }
        }

        return new Resultado(mejorCelulaId, flujoOriginal, menorFlujo);
    }

    private static List<List<Edge>> buildGraph(List<Celula> celulas, double d, Map<Integer, Integer> cellIdToIndex, int source, int sink) {
        int n = celulas.size();
        int totalNodes = n + 2;
        List<List<Edge>> graph = new ArrayList<>(totalNodes);
        for (int i = 0; i < totalNodes; i++) {
            graph.add(new ArrayList<>());
        }

        for (int i = 0; i < n; i++) {
            Celula c = celulas.get(i);
            cellIdToIndex.put(c.id, i);
            indexToCellId.put(i, c.id);
        }

        for (int i = 0; i < n; i++) {
            Celula c1 = celulas.get(i);
            if (c1.tipo == 1) { // Type 1 cells connect to source
                addEdge(graph, source, i, Integer.MAX_VALUE);
            }
            if (c1.tipo == 3) { // Type 3 cells connect to sink
                addEdge(graph, i, sink, Integer.MAX_VALUE);
            }
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    Celula c2 = celulas.get(j);
                    if (puedenConectarsePorTipo(c1, c2) && puedenConectarse(c1, c2, d)) {
                        int peptidosCompartidos = contarPeptidosCompartidos(c1, c2);
                        if (peptidosCompartidos > 0) {
                            addEdge(graph, i, j, peptidosCompartidos);
                        }
                    }
                }
            }
        }
        return graph;
    }

    private static void addEdge(List<List<Edge>> graph, int from, int to, int capacity) {
        Edge forward = new Edge(to, capacity);
        Edge backward = new Edge(from, 0);
        forward.reverse = backward;
        backward.reverse = forward;
        graph.get(from).add(forward);
        graph.get(to).add(backward);
    }

    static private int bloqueadaIndexGlobal = -1;

    private static int maxFlow(List<List<Edge>> originalGraph, int source, int sink, int bloqueadaIndex) {
        int n = originalGraph.size();
        bloqueadaIndexGlobal = bloqueadaIndex;
        // Deep copy the graph to preserve original capacities
        List<List<Edge>> graph = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            graph.add(new ArrayList<>());
            for (Edge e : originalGraph.get(i)) {
                Edge edgeCopy = new Edge(e.to, e.capacity);
                graph.get(i).add(edgeCopy);
            }
        }

        // Reconstruct reverse edges for the copied graph
        for (int i = 0; i < n; i++) {
            for (Edge e : graph.get(i)) {
                for (Edge rev : graph.get(e.to)) {
                    if (rev.to == i) {
                        e.reverse = rev;
                    }
                }
            }
        }

        if (bloqueadaIndex >= 0) {
            // Remove all edges to and from the blocked node
            graph.get(bloqueadaIndex).clear();
            for (List<Edge> edges : graph) {
                edges.removeIf(e -> e.to == bloqueadaIndex);
            }
        }

        int flow = 0;
        int[] level = new int[n];
        while (bfs(graph, source, sink, level)) {
            int[] ptr = new int[n];
            int pushed;
            while ((pushed = dfs(graph, source, sink, ptr, level, Integer.MAX_VALUE)) > 0) {
                flow += pushed;
            }
        }
        return flow;
    }

    private static boolean bfs(List<List<Edge>> graph, int source, int sink, int[] level) {
        Arrays.fill(level, -1);
        Queue<Integer> queue = new LinkedList<>();
        queue.offer(source);
        level[source] = 0;
        while (!queue.isEmpty()) {
            int v = queue.poll();
            for (Edge e : graph.get(v)) {
                if (e.capacity > 0 && level[e.to] < 0) {
                    level[e.to] = level[v] + 1;
                    queue.offer(e.to);
                }
            }
        }
        return level[sink] >= 0;
    }

    private static int dfs(List<List<Edge>> graph, int v, int sink, int[] ptr, int[] level, int flow) {
        if (v == sink || flow == 0) {
            return flow;
        }
        for (; ptr[v] < graph.get(v).size(); ptr[v]++) {
            Edge e = graph.get(v).get(ptr[v]);
            if (level[e.to] == level[v] + 1 && e.capacity > 0) {
                int pushed = dfs(graph, e.to, sink, ptr, level, Math.min(flow, e.capacity));
                if (pushed > 0) {
                    e.capacity -= pushed;
                    e.reverse.capacity += pushed;
                    if (bloqueadaIndexGlobal == -1) {
                        maxFlowVertices.add(indexToCellId.get(v));
                    }
                    return pushed;
                }
            }
        }
        return 0;
    }

    private static class Edge {
        int to;
        int capacity;
        Edge reverse;

        public Edge(int to, int capacity) {
            this.to = to;
            this.capacity = capacity;
        }
    }

    private static boolean puedenConectarsePorTipo(Celula c1, Celula c2) {
        if (c1.tipo == 1)
            return c2.tipo == 2;
        if (c1.tipo == 2)
            return c2.tipo == 2 || c2.tipo == 3;
        return false;
    }

    private static boolean puedenConectarse(Celula c1, Celula c2, double d) {
        double distancia = Math.hypot(c1.posX - c2.posX, c1.posY - c2.posY);
        return distancia <= d;
    }

    private static int contarPeptidosCompartidos(Celula c1, Celula c2) {
        Set<String> peptidos1 = new HashSet<>(c1.peptidos);
        int count = 0;
        for (String p : c2.peptidos) {
            if (peptidos1.contains(p)) {
                count++;
            }
        }
        return count;
    }
}
