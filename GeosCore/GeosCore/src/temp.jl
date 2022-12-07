function a()
	for i = 1:10
		for j = 1:10
			println("< $(i):$(j)")
			
			if j % 2 == 0
				@goto ex
			end
			
			println("> $(i):$(j)")
		end
		@label ex
	end
end

a()
