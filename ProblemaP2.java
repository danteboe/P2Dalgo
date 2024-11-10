//Leer archivo

//Armar mapa de nodos

//Hacer BFS para encontrar nodos alcanzables

//Encontrar capacidad entre neuronas

//Ford-Fulkerson para hallar flujo máximo

//Encontrar célula calculadora con mayor flujo
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;

class ProblemaP2 {
    public static void main(String[] args) throws Exception {
        ProblemaP2 instancia = new ProblemaP2();
        try (
                InputStreamReader is = new InputStreamReader(System.in);
                BufferedReader br = new BufferedReader(is);) {
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

    private static Resultado resolverCaso(List<Celula> celulas, double d) {
        int flujoOriginal = calcularFlujoMaximo(celulas, d, -1);

        int mejorCelulaId = -1;
        int menorFlujo = Integer.MAX_VALUE;

        for (int i = celulas.size() - 1; i >= 0; i--) {
            Celula celula = celulas.get(i);
            if (celula.tipo == 2) { // Calculadora
                int flujoActual = calcularFlujoMaximo(celulas, d, celula.id);
                if (flujoActual <= menorFlujo) {
                    menorFlujo = flujoActual;
                    mejorCelulaId = celula.id;
                }
            }
        }

        return new Resultado(mejorCelulaId, flujoOriginal, menorFlujo);
    }

    private static int calcularFlujoMaximo(List<Celula> celulas, double d, int celulaBloqueada) {
        int n = celulas.size();
        List<List<Integer>> grafo = new ArrayList<>();
        for (int i = 0; i < n + 2; i++) {
            grafo.add(new ArrayList<>());
        }

        int[][] capacidad = new int[n + 2][n + 2];
        int source = n;
        int sink = n + 1;

        for (int i = 0; i < n; i++) {
            Celula c1 = celulas.get(i);
            if (c1.id == celulaBloqueada)
                continue;

            if (c1.tipo == 1) {
                grafo.get(source).add(i);
                grafo.get(i).add(source);
                capacidad[source][i] = Integer.MAX_VALUE;
            }

            if (c1.tipo == 3) {
                grafo.get(i).add(sink);
                grafo.get(sink).add(i);
                capacidad[i][sink] = Integer.MAX_VALUE;
            }

            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;

                Celula c2 = celulas.get(j);
                if (c2.id == celulaBloqueada)
                    continue;

                if (!puedenConectarsePorTipo(c1, c2))
                    continue;

                if (puedenConectarse(c1, c2, d)) {
                    int peptidosCompartidos = contarPeptidosCompartidos(c1, c2);
                    if (peptidosCompartidos > 0) {
                        grafo.get(i).add(j);
                        grafo.get(j).add(i);
                        capacidad[i][j] = peptidosCompartidos;
                        capacidad[j][i] = peptidosCompartidos;
                    }
                }
            }
        }

        return fordFulkerson(grafo, capacidad, source, sink);
    }

    private static boolean puedenConectarsePorTipo(Celula c1, Celula c2) {
        // Iniciadora: 1
        // Calculadora: 2
        // Ejecutora: 3

        // Iniciadora solo puede enviar calculadora
        if (c1.tipo == 1)
            return c2.tipo == 2;

        // Calculadora puede enviar a tipo 2 y 3
        if (c1.tipo == 2)
            return c2.tipo == 2 || c2.tipo == 3;

        // Ejecutor asolo puede recibir, no enviar
        if (c1.tipo == 3)
            return false;

        return false;
    }

    private static boolean puedenConectarse(Celula c1, Celula c2, double d) {
        double distancia = Math.sqrt(Math.pow(c1.posX - c2.posX, 2) + Math.pow(c1.posY - c2.posY, 2));
        return distancia <= d;
    }

    private static int contarPeptidosCompartidos(Celula c1, Celula c2) {
        Set<String> peptidos1 = new HashSet<>(c1.peptidos);
        return (int) c2.peptidos.stream().filter(peptidos1::contains).count();
    }

    private static int fordFulkerson(List<List<Integer>> grafo, int[][] capacidad, int source, int sink) {
        int n = grafo.size();
        int flujoMaximo = 0;
        int[][] flujo = new int[n][n];

        while (true) {
            int[] padre = new int[n];
            Arrays.fill(padre, -1);
            Queue<Integer> cola = new LinkedList<>();
            cola.add(source);
            padre[source] = source;

            while (!cola.isEmpty() && padre[sink] == -1) {
                int u = cola.poll();
                for (int v : grafo.get(u)) {
                    if (padre[v] == -1 && capacidad[u][v] - flujo[u][v] > 0) {
                        padre[v] = u;
                        cola.add(v);
                    }
                }
            }

            if (padre[sink] == -1)
                break;

            int flujoPath = Integer.MAX_VALUE;
            for (int v = sink; v != source; v = padre[v]) {
                int u = padre[v];
                flujoPath = Math.min(flujoPath, capacidad[u][v] - flujo[u][v]);
            }

            for (int v = sink; v != source; v = padre[v]) {
                int u = padre[v];
                flujo[u][v] += flujoPath;
                flujo[v][u] -= flujoPath;
            }

            flujoMaximo += flujoPath;
        }

        return flujoMaximo;
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
}
