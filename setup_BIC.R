library(gtools)
# search recebe char caractere e father inteiro que
# é o endereço do nó que terá um nó filho rotulado com char.
# Ele procura entre os filhos de nó father um nó rotulado
# com char.
# Se encontrar ele retorna uma lista de inteiros com os endereços
# do nó irmão anterior, o próprio nó e o nó irmão posterior repectivamente.
# Se não ele retorna a lista com três NA`s
search <- defmacro(tree, char, father, expr={
	prev_sibling <- NA
	child <- NA
	next_sibling <- NA
	if(!is.na(tree$mat[father, 3])){
		next_sibling <- as.numeric(tree$mat[father, 3])
		repeat{
			prev_sibling <- child
			child <- next_sibling
			next_sibling <- as.numeric(tree$mat[next_sibling, 2])
			if(tree$mat[child, 4] >= char | is.na(next_sibling)) break()
		}
		if(tree$mat[child, 4] != char){
			if(tree$mat[child, 4] > char) next_sibling <- child
			else prev_sibling <- child
			child <- NA
		}
	}
})

# insert recebe char caractere e father inteiro que 
# é o endereço do nó que terá um nó filho rotulado com char.
# Se não existir um nó rotulado com char entre o filhos do father
# um novo nó é criado.
new_node <- defmacro(tree, expr={
	if(tree$next_free_node == tree$total_row + 1){
		tree$total_row <- tree$total_row + 1024
		tree$mat <- rbind(tree$mat, matrix(nrow = 1024, ncol = 8))
	}
	node <- tree$next_free_node
	tree$next_free_node <- tree$next_free_node + 1
})
insert <- defmacro(tree, char, father, expr={
	search(tree, char, father)
	if(is.na(child)){
		new_node(tree)
		tree$mat[node, 1] <- father
		tree$mat[node, 2] <- next_sibling
		tree$mat[node, 4] <- char
		tree$mat[node, 7] <- 0
		if(is.na(prev_sibling)) tree$mat[father, 3] <- node
		else tree$mat[prev_sibling, 2] <- node
		father <- node
	}
	else father <- child
})

degree <- function(pre_node, tree){
	s <- 0
	node <- as.numeric(tree$mat[pre_node, 3])
	while(!is.na(node)){
		s <- s + 1
		node <- as.numeric(tree$mat[node, 2])
	}
	return(s)
}

calculate_probabilities <- function(pre_node, tree){
	if(is.na(pre_node)) return(tree)
	sibling <- pre_node
	occ <- 0
	while(!is.na(sibling)){
		occ <- as.numeric(tree$mat[sibling, 7]) + occ
		sibling <- as.numeric(tree$mat[sibling, 2])
	}
	while(!is.na(pre_node)){
		# set_degrees_freedom
		tree$mat[pre_node, 6] <- 0
		s <- 0
		node <- as.numeric(tree$mat[pre_node, 3])
		while(!is.na(node)){
			s <- s + 1
			node <- as.numeric(tree$mat[node, 2])
		}
		tree$mat[pre_node, 6] <- s
		# set_probability
		tree$mat[pre_node, 5] <- as.numeric(tree$mat[pre_node, 7]) / occ
		# Vai para o filho
		tree <- calculate_probabilities(as.numeric(tree$mat[pre_node, 3]), tree)
		# Vai para o irmão
		pre_node <- as.numeric(tree$mat[pre_node, 2])
	}
	return(tree)
}

calculate_ell <- function(pre_node, tree){
	if(is.na(pre_node)) return(tree)
	value <- 0
	node <- as.numeric(tree$mat[pre_node, 3])
	while(!is.na(node)){
		value <- value + as.numeric(tree$mat[node, 7]) * log(as.numeric(tree$mat[node, 5]))
		node <- as.numeric(tree$mat[node, 2])
	}
	tree$mat[pre_node, 8] <- value
	tree <- calculate_ell(as.numeric(tree$mat[pre_node, 3]), tree)
	tree <- calculate_ell(as.numeric(tree$mat[pre_node, 2]), tree)
	return(tree)
}

get_pre_node <- function(suf_node, tree){
	pre_node <- 1
	while(!is.na(tree$mat[suf_node, 1])){
		pre_node <- as.numeric(tree$mat[pre_node, 3])
		while(tree$mat[pre_node, 4] != tree$mat[suf_node, 4]) pre_node <- as.numeric(tree$mat[pre_node, 2])
		suf_node <- as.numeric(tree$mat[suf_node, 1])
	}
	return(pre_node)
}

set_mate <- function(tree){
	last_suf_node <- tree$next_free_node - 1
	for(suf_node in tree$suf : last_suf_node) tree$mat[suf_node, 8] <- get_pre_node(suf_node, tree)
	return(tree)
}

# setup_BIC é um função que recebe sample, uma cadeia de caracteres,
# e depth, número de nós de profundidade na árvore de sufixo.
# Ela cria tree, uma matrix que representa pref e suf duas árvores,
# que serão usadas para o BIC, fazendo as configuração certas para tal fim.
setup_BIC <- function(sample, depth){
	# Criando árvores
	tree = list(mat = matrix(nrow = 1024, ncol = 8), next_free_node = 1, total_row = 1024, suf = NA)
	# Criando pre tree
	tree$mat[tree$next_free_node, 6] <- 0
	tree$mat[tree$next_free_node, 7] <- size_sample <- nchar(sample)
	tree$next_free_node <- tree$next_free_node + 1
	i <- 1
	while(i <= size_sample){
		father <- 1
		for(char in unlist(strsplit(substr(sample, i, i + depth), ""))){
			insert(tree, char, father)
			tree$mat[father, 7] <- as.numeric(tree$mat[father, 7]) + 1
		}
		i <- i + 1
	}
	tree <- calculate_probabilities(as.numeric(tree$mat[1, 3]), tree)
	tree <- calculate_ell(1, tree)
	# Criando suf tree
	new_node(tree)
	tree$suf <- node
	i <- 1
	depth <- depth - 1
	while(i <= size_sample){
		father <- tree$suf
		for(char in rev(unlist(strsplit(substr(sample, i - depth, i), "")))) insert(tree, char, father)
		i <- i + 1
	}
	tree <- set_mate(tree)
	return(tree)
}