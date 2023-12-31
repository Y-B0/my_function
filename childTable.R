childTable <- function(x, vars = NULL, opts = NULL, ...) {
  names_x <- names(x)
  if (is.null(vars)) stop("'vars' must be specified!")
  pos <- match(vars, names_x)
  pos <- pos[pos <= ncol(x)] + 1
  rownames(x) <- NULL
  if (nrow(x) > 0) x <- cbind(' ' = '&oplus;', x)
  # options
  opts <- c(
    opts,
    scrollX = TRUE,
    list(
      columnDefs = list(
        list(visible = FALSE, targets = c(0, pos)),
        list(orderable = FALSE, className = 'details-control', targets = 1),
        list(className = 'dt-left', targets = 1:3),
        list(className = 'dt-right', targets = 4:ncol(x))
      )
    ))
  datatable(
    x,
    ...,
    escape = F,
    options = opts,
    callback = JS(.callback2(x = x, pos = c(0, pos)))
  )
}
.callback2 <- function(x, pos = NULL) {
  part1 <- "table.column(1).nodes().to$().css({cursor: 'pointer'});"
  part2 <- .child_row_table2(x, pos = pos)
  part3 <-
    "
  table.on('click', 'td.details-control', function() {
  var td = $(this), row = table.row(td.closest('tr'));
  if (row.child.isShown()) {
  row.child.hide();
  td.html('&oplus;');
  } else {
  row.child(format(row.data())).show();
  td.html('&ominus;');
  }
  });"
  
  paste(part1, part2, part3)
}
.child_row_table2 <- function(x, pos = NULL) {
  names_x <- paste0(names(x), ":")
  text <- "
  var format = function(d) {
  text = '<div><table>' +
  "
  
  for (i in seq_along(pos)) {
    text <- paste(text, glue::glue(
      "'<tr>' +
      '<td>' + '{names_x[pos[i]]}' + '</td>' +
      '<td style="background-color:#eee; word-wrap:break-word;word-break:break-all; ">' + d[{pos[i]}] + '</td>' +
        '</tr>' + " ))
  }
  paste0(text,
         "'</table></div>'
      return text;};"
  )
}