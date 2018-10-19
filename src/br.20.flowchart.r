require(DiagrammeR) 

dirp = '~/briggs/data'
setwd(dirp)
getwd()

graph = grViz("
digraph nicegraph {

  # graph, node, and edge definitions
  graph [compound = true, nodesep = .5, ranksep = .25,
         color = crimson]

  node [fontname = Helvetica, fontcolor = darkslategray,
        shape = rectangle, fixedsize = true, width = 1,
        color = darkslategray]

  edge [color = grey, arrowhead = none, arrowtail = none]

  # subgraph for R information
  subgraph cluster0 {
    node [fixedsize = true, width = 3]
    '@@1-1' -> '@@1-2' -> '@@1-3' -> '@@1-4'
    '@@1-4' -> '@@1-5' -> '@@1-6' -> '@@1-7'
  }

  # subgraph for RStudio information
  subgraph cluster1 {
    node [fixedsize = true, width = 3]
    '@@2' -> '@@3'
  }

  Information             [width = 1.5]
  Information -> R
  Information -> RStudio
  R -> '@@1-1'            [lhead = cluster0]
  RStudio -> '@@2'        [lhead = cluster1]

}

[1]: paste0(names(R.Version())[1:7], ':\\n ', R.Version()[1:7])
[2]: paste0('RStudio version:\\n ', '5.6')
[3]: paste0('Current program mode:\\n ', '6.9')

")

fp = file.path(dirp, 'flowchart.pdf')
graph %>%
    export_graph(file_name = fp, file_type = "PDF", width = 8, height = 6)

pdf(fp, width = 8, height = 6)
print(graph)
dev.off()
